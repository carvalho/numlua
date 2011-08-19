/* {=================================================================
 *
 * stat.c
 * Statistical routines for NumericLua
 * Luis Carvalho (lexcarvalho@gmail.com)
 * See Copyright Notice in numlua.h
 *
 * ==================================================================} */

#include <lua.h>
#include <lauxlib.h>
#include <float.h>
#include <math.h>
#include "numlua.h"
#include "cdflib.h"

#define MAXITER (1000) /* used in dchisq */
#define LBOUND (1.0e-20) /* ditto */
#define SQRT2PI (2.506628274631) /* sqrt(2*pi) */
/* From Catherine Loader's library (dbinom) */
#define FORCE_INT(x) floor((x) + 0.5)
/* used in stirlerr: */
#define HF_LG_PIx2 0.918938533204672741780329736406 /* 0.5*log(2*pi) */
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508 /* 1/1188 */
/* error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0. */
static double sferr_halves[31] = {
  0.0, /* n=0 - wrong, place holder only */
  0.1534264097200273452913848,  /* 0.5 */
  0.0810614667953272582196702,  /* 1.0 */
  0.0548141210519176538961390,  /* 1.5 */
  0.0413406959554092940938221,  /* 2.0 */
  0.03316287351993628748511048, /* 2.5 */
  0.02767792568499833914878929, /* 3.0 */
  0.02374616365629749597132920, /* 3.5 */
  0.02079067210376509311152277, /* 4.0 */
  0.01848845053267318523077934, /* 4.5 */
  0.01664469118982119216319487, /* 5.0 */
  0.01513497322191737887351255, /* 5.5 */
  0.01387612882307074799874573, /* 6.0 */
  0.01281046524292022692424986, /* 6.5 */
  0.01189670994589177009505572, /* 7.0 */
  0.01110455975820691732662991, /* 7.5 */
  0.010411265261972096497478567, /* 8.0 */
  0.009799416126158803298389475, /* 8.5 */
  0.009255462182712732917728637, /* 9.0 */
  0.008768700134139385462952823, /* 9.5 */
  0.008330563433362871256469318, /* 10.0 */
  0.007934114564314020547248100, /* 10.5 */
  0.007573675487951840794972024, /* 11.0 */
  0.007244554301320383179543912, /* 11.5 */
  0.006942840107209529865664152, /* 12.0 */
  0.006665247032707682442354394, /* 12.5 */
  0.006408994188004207068439631, /* 13.0 */
  0.006171712263039457647532867, /* 13.5 */
  0.005951370112758847735624416, /* 14.0 */
  0.005746216513010115682023589, /* 14.5 */
  0.005554733551962801371038690  /* 15.0 */
};



/* {=====================================================================
 *    Factors
 * ======================================================================} */

typedef struct {
  int size;
  int nlevels;
  unsigned char map[1];
} nl_Factor;

#define MAXLEVELS 255
#define LEVELS 4
#define MAP 5
static int stat_factor (lua_State *L) {
  int i, n = (int) lua_rawlen(L, 1);
  nl_Factor *f;
  unsigned char m, l = 0; /* # levels */
  lua_settop(L, 1);
  luaL_argcheck(L, n > 0, 1, "length must be positive");
  f = (nl_Factor *) lua_newuserdata(L, sizeof(nl_Factor) + n - 1);
  f->size = n;
  f->nlevels = 0;
  lua_pushvalue(L, -1);
  lua_newtable(L); /* levels */
  lua_newtable(L); /* map */
  for (i = 0; i < n; i++) {
    lua_pushinteger(L, i + 1);
    lua_gettable(L, 1); /* t[i] */
    lua_pushvalue(L, -1);
    lua_gettable(L, MAP); /* map[t[i]] */
    if (lua_isnil(L, -1)) { /* add level? */
      if (f->nlevels == MAXLEVELS)
        luaL_error(L, "maximum number of levels exceeded");
      m = (unsigned char) ++f->nlevels;
      lua_pop(L, 1); /* nil */
      lua_pushvalue(L, -1); /* t[i] */
      lua_pushinteger(L, (int) m);
      lua_settable(L, MAP); /* map[t[i]] = m */
      lua_rawseti(L, LEVELS, (int) m); /* levels[m] = t[i] */
    }
    else {
      m = (unsigned char) lua_tointeger(L, -1);
      lua_pop(L, 2); /* t[i], map[t[i]] */
    }
    f->map[i] = m - 1; /* zero based */
  }
  lua_pop(L, 1); /* map */
  lua_rawset(L, lua_upvalueindex(1)); /* metatable[f] = levels */
  lua_pushvalue(L, lua_upvalueindex(1)); /* metatable */
  lua_setmetatable(L, 2);
  return 1;
}


static int factor__tostring (lua_State *L) {
  lua_pushfstring(L, "factor: %p", lua_touserdata(L, 1));
  return 1;
}

static int factor__len (lua_State *L) {
  nl_Factor *f = (nl_Factor *) lua_touserdata(L, 1);
  nl_Matrix *l = nl_pushmatrix(L, 0, 1, &f->nlevels, 1, f->nlevels, NULL);
  int i;
  for (i = 0; i < l->size; i++) l->data[i] = 0;
  for (i = 0; i < f->size; i++) l->data[f->map[i]]++;
  return 1;
}

static int factor__call (lua_State *L) {
  nl_Factor *f = (nl_Factor *) lua_touserdata(L, 1);
  if (lua_gettop(L) == 1) /* levels? */
    lua_rawget(L, lua_upvalueindex(1));
  else {
    int k = lua_tointeger(L, 2);
    if (k < 1 || k > f->size)
      lua_pushnil(L);
    else
      lua_pushinteger(L, (int) f->map[k - 1] + 1);
  }
  return 1;
}

static int factor__index (lua_State *L) {
  nl_Factor *f = (nl_Factor *) lua_touserdata(L, 1);
  if (lua_isnumber(L, 2)) {
    int k = lua_tointeger(L, 2);
    if (k < 1 || k > f->size)
      lua_pushnil(L);
    else {
      lua_settop(L, 1);
      lua_rawget(L, lua_upvalueindex(1));
      lua_rawgeti(L, -1, (int) f->map[k - 1] + 1);
    }
  }
  else /* meta lookup? */
    lua_rawget(L, lua_upvalueindex(2));
  return 1;
}


static int factor_fold (lua_State *L) {
  nl_Factor *f = (nl_Factor *) lua_touserdata(L, 1);
  nl_Matrix *r, *m = nl_checkmatrix(L, 2);
  nl_Complex c = nl_optcomplex(L, 4, 0);
  int i;
  luaL_argcheck(L, m->size == f->size, 2, "inconsistent sizes");
  luaL_argcheck(L, !m->section, 2, "sections are not allowed");
  luaL_argcheck(L, lua_type(L, 3) == LUA_TFUNCTION, 3,
      "function expected");
  lua_settop(L, 4);
  /* result vector: */
  r = nl_pushmatrix(L, m->iscomplex, 1, &f->nlevels, 1, f->nlevels, NULL);
  /* fold: */
  if (r->iscomplex) {
    for (i = 0; i < f->nlevels; i++) /* init result */
      CPX(r->data)[i] = c;
    for (i = 0; i < f->size; i++) {
      lua_pushvalue(L, 3); /* function */
      nl_pushcomplex(L, CPX(r->data)[f->map[i]]);
      nl_pushcomplex(L, CPX(m->data)[i]);
      lua_call(L, 2, 1);
      CPX(r->data)[f->map[i]] = nl_optcomplex(L, -1, 0);
      lua_pop(L, 1);
    }
  }
  else {
    for (i = 0; i < f->nlevels; i++) /* init result */
      r->data[i] = creal(c);
    for (i = 0; i < f->size; i++) {
      lua_pushvalue(L, 3); /* function */
      lua_pushnumber(L, r->data[f->map[i]]);
      lua_pushnumber(L, m->data[i]);
      lua_call(L, 2, 1);
      r->data[f->map[i]] = luaL_optnumber(L, -1, 0);
      lua_pop(L, 1);
    }
  }
  return 1;
}

static int factor_partition (lua_State *L) {
  nl_Factor *f = (nl_Factor *) lua_touserdata(L, 1);
  nl_Matrix *m = nl_checkmatrix(L, 2);
  nl_Matrix *r[MAXLEVELS];
  int i, l[MAXLEVELS];
  luaL_argcheck(L, m->size == f->size, 2, "inconsistent sizes");
  luaL_argcheck(L, !m->section, 2, "sections are not allowed");
  /* process level counts */
  for (i = 0; i < f->nlevels; i++) l[i] = 0;
  for (i = 0; i < f->size; i++) l[f->map[i]]++;
  /* create partition vectors */
  lua_newtable(L);
  for (i = 0; i < f->nlevels; i++) { /* init vectors and counts */
    r[i] = nl_pushmatrix(L, m->iscomplex, 1, &l[i], 1, l[i], NULL);
    l[i] = 0;
    lua_rawseti(L, -2, i + 1);
  }
  if (m->iscomplex) {
    for (i = 0; i < f->size; i++)
      CPX(r[f->map[i]]->data)[l[i]++] = CPX(m->data)[i];
  }
  else {
    for (i = 0; i < f->size; i++) {
      int k = f->map[i];
      r[k]->data[l[k]++] = m->data[i];
    }
  }
  return 1;
}

static int factor_design (lua_State *L) {
  nl_Factor *f = (nl_Factor *) lua_touserdata(L, 1);
  int i, j, ref = luaL_optinteger(L, 2, 0);
  nl_Matrix *m;
  lua_Number *e;
  luaL_argcheck(L, ref >= 0 || ref <= f->nlevels, 2,
      "invalid reference class");
  m = nl_pushmatrix(L, 0, 2, NULL, 1, f->nlevels * f->size, NULL);
  m->dim[0] = f->size; m->dim[1] = f->nlevels;
  ref--; /* zero based */
  e = m->data;
  for (i = 0; i < f->nlevels; i++) {
    if (i == ref) {
      for (j = 0; j < f->size; j++) *e++ = 1; /* intercept */
    }
    else {
      for (j = 0; j < f->size; j++)
        *e++ = (i == f->map[j]) ? 1 : -(ref == f->map[j]);
    }
  }
  return 1;
}


/* {=====================================================================
 *    Auxiliary
 * ======================================================================} */

/* Failsafe: execution shouldn't reach here (!), since most errors are checked
 * out by specific check_xxx routines; the only expected error is when status
 * == 10 */
static void check_status (lua_State *L, int status, lua_Number bound) {
  if (status == 1)
    luaL_error(L, "result lower than search bound: %f", bound);
  if (status == 2)
    luaL_error(L, "result higher than search bound: %f", bound);
  if (status < 0)
    luaL_error(L, "out of range on parameter %d: %f", -status, bound);
  if (status == 10)
    luaL_error(L, "error in cumgam: %d", status);
}


/* {=======   Catherine Loader's dbinom   =======}
 *
 * The code below -- bd0, stirlerr, and dbinom_raw -- is adapted from
 * Catherine Loader <c at herine dot net>
 * http://www.herine.net/stat/software/dbinom.html
 *
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, with the exceptions noted below,
 * and provided that this entire notice is included in all copies of any
 * software which is or includes a copy or modification of this software
 * and in all copies of the supporting documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR NOR LUCENT TECHNOLOGIES
 * MAKE ANY REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE
 * MERCHANTABILITY OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/* bd0: computes bd0(x, np) = x log(x/np) + np - x */
static lua_Number bd0 (lua_Number x, lua_Number np) {
  if (fabs(x - np) < 0.1 * (x + np)) {
    lua_Number ej, s, s1, v;
    int j;
    s = (x - np);
    s *= v = s / (x + np);
    ej = 2 * x * v;
    v *= v;
    for (j = 1; ; j++) {
      ej *= v;
      s1 = s + ej / ((j << 1) + 1);
      if (s1 == s) return s1;
      s = s1;
    }
  }
  return x * log(x / np) + np - x;
}

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n ) */
static lua_Number stirlerr (lua_Number n) {
  lua_Number nn;
  if (n < 15.0) {
    nn = 2.0 * n;
    if (nn == FORCE_INT(nn)) return sferr_halves[(int) nn];
    return lgamma(n + 1.0) - (n + 0.5) * log(n) + n - HF_LG_PIx2;
  }
  nn = n * n;
  if (n > 500) return (S0 - S1 / nn) / n;
  if (n > 80) return (S0 - (S1 - S2 / nn) / nn) / n;
  if (n > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

static lua_Number dbinom_raw (lua_Number x, lua_Number n, lua_Number p,
    lua_Number q) {
  lua_Number f, lc;
  if (p == 0) return (x == 0) ? 1 : 0;
  if (q == 0) return (x == n) ? 1 : 0;
  if (x == 0)
    return exp((p < 0.1) ? -bd0(n, n * q) - n * p : n * log(q));
  if (x == n)
    return exp((q < 0.1) ? -bd0(n, n * p) - n * q : n * log(p));
  if ((x < 0) || (x > n)) return 0;
  lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x)
    - bd0(x, n * p) - bd0(n - x, n * q);
  f = (2 * M_PI * x * (n - x)) / n;
  return exp(lc) / sqrt(f);
}

static lua_Number dhyper_raw (lua_Number x, lua_Number r, lua_Number b,
    lua_Number n) {
  lua_Number p, q, p1, p2, p3;
  if (x < 0) return 0;
  x = FORCE_INT(x);
  r = FORCE_INT(r);
  b = FORCE_INT(b);
  n = FORCE_INT(n);
  if (n == 0) return (x == 0) ? 1 : 0;
  p = n / (r + b);
  q = (r + b - n) / (r + b);
  p1 = dbinom_raw(x, r, p, q);
  p2 = dbinom_raw(n - x, b, p, q);
  p3 = dbinom_raw(n, r + b, p, q);
  return p1 * p2 / p3;
}

/* Code adapted from R (2.4.0) */
static lua_Number pdhyper (lua_Number x, lua_Number r, lua_Number b,
    lua_Number n) {
  lua_Number s = 0, t = 1;
  while (x > 0 && t >= DBL_EPSILON * s) {
    t *= x * (b - n + x) / (n + 1 - x) / (r + 1 - x);
    s += t;
    x--;
  }
  return 1 + s;
}


/* {=====================================================================
 *    Functions
 * ======================================================================} */

/* {=======   Beta   =======} */

static void check_beta (lua_State *L, int which, lua_Number x,
    lua_Number a, lua_Number b) {
  (void) which;
  luaL_argcheck(L, x >= 0 && x <= 1, 1, "out of range");
  luaL_argcheck(L, a >= 0, 2, "non-negative value expected");
  luaL_argcheck(L, b >= 0, 3, "non-negative value expected");
}

static int stat_dbeta (lua_State *L) {
  /* stack should contain x, a and b */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number a = luaL_checknumber(L, 2);
  lua_Number b = luaL_checknumber(L, 3);
  check_beta(L, 1, x, a, b);
  lua_pushnumber(L, (x == 0 || x == 1) ? 0 :
      exp((a - 1) * log(x) + (b - 1) * log(1 - x) - dlnbet(&a, &b)));
  return 1;
}

static int stat_pbeta (lua_State *L) {
  /* stack should contain x, a and b */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number a = luaL_checknumber(L, 2);
  lua_Number b = luaL_checknumber(L, 3);
  lua_Number p, q, y, bound;
  int which = 1;
  int status;
  check_beta(L, 1, x, a, b);
  y = 1 - x;
  cdfbet(&which, &p, &q, &x, &y, &a, &b, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qbeta (lua_State *L) {
  /* stack should contain x, a and b */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number a = luaL_checknumber(L, 2);
  lua_Number b = luaL_checknumber(L, 3);
  lua_Number x;
  check_beta(L, 2, p, a, b);
  if (p == 0 || p == 1) x = p;
  else {
    lua_Number q, y, bound;
    int which = 2;
    int status;
    q = 1 - p;
    cdfbet(&which, &p, &q, &x, &y, &a, &b, &status, &bound);
    check_status(L, status, bound);
  }
  lua_pushnumber(L, x);
  return 1;
}


/* {=======   Binomial   =======} */

static void check_binom (lua_State *L, int which, lua_Number x,
    lua_Number xn, lua_Number pr) {
  luaL_argcheck(L, ((which == 1 && (x >= 0 && x <= xn)) /* x */
      || (which == 2 && (x >= 0 && x <= 1))), /* p */
      1, "out of range");
  luaL_argcheck(L, xn >= 0, 2, "non-negative value expected");
  luaL_argcheck(L, pr >= 0 && pr <= 1, 3, "out of range");
}

static int stat_dbinom (lua_State *L) {
  /* stack should contain s, xn, pr */
  lua_Number s = luaL_checknumber(L, 1);
  lua_Number xn = luaL_checknumber(L, 2);
  lua_Number pr = luaL_checknumber(L, 3);
  check_binom(L, 1, s, xn, pr);
  s = FORCE_INT(s);
  xn = FORCE_INT(xn);
  lua_pushnumber(L, dbinom_raw(s, xn, pr, 1 - pr));
  return 1;
}

static int stat_pbinom (lua_State *L) {
  /* stack should contain s, xn, pr */
  lua_Number s = luaL_checknumber(L, 1);
  lua_Number xn = luaL_checknumber(L, 2);
  lua_Number pr = luaL_checknumber(L, 3);
  lua_Number p, q, ompr, bound;
  int which = 1;
  int status;
  check_binom(L, 1, s, xn, pr);
  ompr = 1 - pr;
  cdfbin(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qbinom (lua_State *L) {
  /* stack should contain p, xn, pr */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number xn = luaL_checknumber(L, 2);
  lua_Number pr = luaL_checknumber(L, 3);
  lua_Number s;
  int si;
  check_binom(L, 2, p, xn, pr);
  if (p == 0 || p == 1) s = p*xn;
  else {
    lua_Number q = 1 - p;
    lua_Number ompr = 1 - pr;
    lua_Number bound;
    int which = 2;
    int status;
    cdfbin(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
    check_status(L, status, bound);
  }
  lua_number2int(si, s);
  lua_pushinteger(L, si);
  return 1;
}


/* {=======   Chi-square   =======} */

static void check_chisq (lua_State *L, int which, lua_Number x,
    lua_Number df, lua_Number pnonc) {
  luaL_argcheck(L, ((which == 1 && x >= 0)  /* x */
      || (which == 2 && (x >= 0 && x <= 1))), /* p */
      1, "out of range");
  if (pnonc == 0)
    luaL_argcheck(L, df > 0, 2, "positive value expected");
  else
    luaL_argcheck(L, df >= 0, 2, "non-negative value expected");
}

static int stat_dchisq (lua_State *L) {
  /* stack should contain x, df and opt. pnonc */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number df = luaL_checknumber(L, 2);
  lua_Number pnonc = luaL_optnumber(L, 3, 0);
  lua_Number d;
  check_chisq(L, 1, x, df, pnonc);
  /* compute central dchisq */
  d = df / 2;
  d = exp((d - 1) * log(x) - x / 2 - d * M_LN2 - dlngam(&d));
  /* compute non-central if that's the case */
  if (pnonc != 0) { /* non-central? */
    /* evaluate weighted series */
    int i;
    lua_Number t = d *= exp(-pnonc / 2); /* first term */
    for (i = 1; i < MAXITER && d > LBOUND && t > DBL_EPSILON * d; i++)
      d += t *= x * pnonc / (2 * i * (df + 2 * (i - 1)));
  }
  lua_pushnumber(L, d);
  return 1;
}

static int stat_pchisq (lua_State *L) {
  /* stack should contain x, df and opt. pnonc */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number df = luaL_checknumber(L, 2);
  lua_Number pnonc = luaL_optnumber(L, 3, 0);
  lua_Number p, q, bound;
  int which = 1;
  int status;
  check_chisq(L, 1, x, df, pnonc);
  if (pnonc == 0) /* central? */
    cdfchi(&which, &p, &q, &x, &df, &status, &bound);
  else /* non-central */
    cdfchn(&which, &p, &q, &x, &df, &pnonc, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qchisq (lua_State *L) {
  /* stack should contain p, df and opt. pnonc */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number df = luaL_checknumber(L, 2);
  lua_Number pnonc = luaL_optnumber(L, 3, 0);
  lua_Number x;
  check_chisq(L, 2, p, df, pnonc);
  if (p == 0 || p == 1) x = (p == 0) ? 0 : HUGE_VAL;
  else {
    lua_Number q = 1 - p;
    lua_Number bound;
    int which = 2;
    int status;
    if (pnonc == 0) /* central? */
      cdfchi(&which, &p, &q, &x, &df, &status, &bound);
    else /* non-central */
      cdfchn(&which, &p, &q, &x, &df, &pnonc, &status, &bound);
    check_status(L, status, bound);
  }
  lua_pushnumber(L, x);
  return 1;
}


/* {=======   Exponential   =======} */

static int stat_dexp (lua_State *L) {
  /* stack should contain x and rate */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number l = luaL_optnumber(L, 2, 1);
  lua_pushnumber(L, exp(-l * x) * l);
  return 1;
}

static int stat_pexp (lua_State *L) {
  /* stack should contain x and rate */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number l = luaL_optnumber(L, 2, 1);
  lua_pushnumber(L, 1 - exp(-l * x));
  return 1;
}

static int stat_qexp (lua_State *L) {
  /* stack should contain p and rate */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number l = luaL_optnumber(L, 2, 1);
  luaL_argcheck(L, p >= 0 && p <= 1, 1, "out of range");
  lua_pushnumber(L, (p < 1) ? -log(1 - p) / l : HUGE_VAL);
  return 1;
}


/* {=======   F   =======} */

static void check_f (lua_State *L, int which, lua_Number x,
    lua_Number dfn, lua_Number dfd) {
  luaL_argcheck(L, ((which == 1 && x >= 0)  /* x */
      || (which == 2 && (x >= 0 && x <= 1))), /* p */
      1, "out of range");
  luaL_argcheck(L, dfn >= 0, 2, "non-negative value expected");
  luaL_argcheck(L, dfd >= 0, 3, "non-negative value expected");
}

static int stat_df (lua_State *L) {
  /* stack should contain f, dfn, dfd */
  lua_Number f = luaL_checknumber(L, 1);
  lua_Number dfn = luaL_checknumber(L, 2);
  lua_Number dfd = luaL_checknumber(L, 3);
  lua_Number df1, df2, r, d;
  check_f(L, 1, f, dfn, dfd);
  df1 = dfn / 2;
  df2 = dfd / 2;
  r = dfn / dfd;
  d = df1 * log(r) + (df1 - 1) * log(f);
  d -= (df1 + df2) * log(1 + r * f);
  d -= dlnbet(&df1, &df2);
  lua_pushnumber(L, exp(d));
  return 1;
}

static int stat_pf (lua_State *L) {
  /* stack should contain f, dfn, dfd and opt. phonc */
  lua_Number f = luaL_checknumber(L, 1);
  lua_Number dfn = luaL_checknumber(L, 2);
  lua_Number dfd = luaL_checknumber(L, 3);
  lua_Number phonc = luaL_optnumber(L, 4, 0);
  lua_Number p, q, bound;
  int which = 1;
  int status;
  check_f(L, 1, f, dfn, dfd);
  if (phonc == 0) /* central? */
    cdff(&which, &p, &q, &f, &dfn, &dfd, &status, &bound);
  else /* non-central */
    cdffnc(&which, &p, &q, &f, &dfn, &dfd, &phonc, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qf (lua_State *L) {
  /* stack should contain p, dfn, dfd and opt. phonc */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number dfn = luaL_checknumber(L, 2);
  lua_Number dfd = luaL_checknumber(L, 3);
  lua_Number phonc = luaL_optnumber(L, 4, 0);
  lua_Number f;
  check_f(L, 2, p, dfn, dfd);
  if (p == 0 || p == 1) f = (p == 0) ? 0 : HUGE_VAL;
  else {
    lua_Number q = 1 - p;
    lua_Number bound;
    int which = 2;
    int status;
    if (phonc == 0) /* central? */
      cdff(&which, &p, &q, &f, &dfn, &dfd, &status, &bound);
    else /* non-central */
      cdffnc(&which, &p, &q, &f, &dfn, &dfd, &phonc, &status, &bound);
    check_status(L, status, bound);
  }
  lua_pushnumber(L, f);
  return 1;
}


/* {=======   Gamma   =======} */

static void check_gamma (lua_State *L, int which, lua_Number x,
    lua_Number shape, lua_Number scale) {
  luaL_argcheck(L, ((which == 1 && x >= 0)  /* x */
      || (which == 2 && (x >= 0 && x <= 1))), /* p */
      1, "out of range");
  luaL_argcheck(L, shape >= 0, 2, "non-negative value expected");
  luaL_argcheck(L, scale >= 0, 3, "non-negative value expected");
}

/* scale here is 1/rate */
static int stat_dgamma (lua_State *L) {
  /* stack should contain x, shape and opt. scale */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number shape = luaL_checknumber(L, 2);
  lua_Number scale = luaL_optnumber(L, 3, 1);
  lua_Number d;
  check_gamma(L, 1, x, shape, scale);
  d = x * scale;
  d = exp(shape * log(d) - d - dlngam(&shape)) / x;
  lua_pushnumber(L, d);
  return 1;
}

static int stat_pgamma (lua_State *L) {
  /* stack should contain x, shape and opt. scale */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number shape = luaL_checknumber(L, 2);
  lua_Number scale = luaL_optnumber(L, 3, 1);
  lua_Number p, q, bound;
  int which = 1;
  int status;
  check_gamma(L, 1, x, shape, scale);
  cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qgamma (lua_State *L) {
  /* stack should contain p, shape and opt. scale */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number shape = luaL_checknumber(L, 2);
  lua_Number scale = luaL_optnumber(L, 3, 1);
  lua_Number x;
  check_gamma(L, 2, p, shape, scale);
  if (p == 0 || p == 1) x = (p == 0) ? 0 : HUGE_VAL;
  else {
    lua_Number q = 1 - p;
    lua_Number bound;
    int which = 2;
    int status;
    cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
    check_status(L, status, bound);
  }
  lua_pushnumber(L, x);
  return 1;
}


/* {=======   Hypergeometric   =======} */

static void check_hyper (lua_State *L, lua_Number x, lua_Number r,
    lua_Number b, lua_Number n) {
  luaL_argcheck(L, x >= 0 && x <= n, 1, "out of range");
  luaL_argcheck(L, r >= 0, 2, "non-negative value expected");
  luaL_argcheck(L, b >= 0, 3, "non-negative value expected");
  luaL_argcheck(L, n >= 0 && n <= b + r, 4, "out of range");
}

static int stat_dhyper (lua_State *L) {
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number r = luaL_checknumber(L, 2);
  lua_Number b = luaL_checknumber(L, 3);
  lua_Number n = luaL_checknumber(L, 4);
  check_hyper(L, x, r, b, n);
  lua_pushnumber(L, dhyper_raw(x, r, b, n));
  return 1;
}

static int stat_phyper (lua_State *L) {
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number r = luaL_checknumber(L, 2);
  lua_Number b = luaL_checknumber(L, 3);
  lua_Number n = luaL_checknumber(L, 4);
  lua_Number d;
  int lower_tail = 1;
  x = FORCE_INT(x);
  r = FORCE_INT(r);
  b = FORCE_INT(b);
  n = FORCE_INT(n);
  check_hyper(L, x, r, b, n);
  if (x * (r + b) > n * r) { /* swap tails? */
    lua_Number s = b;
    b = r;
    r = s;
    x = n - x - 1;
    lower_tail = 0;
  }
  if (x < 0) return 0;
  d = dhyper_raw(x, r, b, n);
  d *= pdhyper(x, r, b, n);
  lua_pushnumber(L, (lower_tail ? d : (1.0 - d)));
  return 1;
}

/* TODO: qhyper */


/* {=======   Negative binomial   =======} */

static void check_nbinom (lua_State *L, int which, lua_Number x,
    lua_Number xn, lua_Number pr) {
  luaL_argcheck(L, ((which == 1 && x >= 0)  /* x */
      || (which == 2 && (x >= 0 && x <= 1))), /* p */
      1, "out of range");
  luaL_argcheck(L, xn >= 0, 2, "non-negative value expected");
  luaL_argcheck(L, pr >= 0 && pr <= 1, 3, "out of range");
}

static int stat_dnbinom (lua_State *L) {
  /* stack should contain s, xn, pr */
  lua_Number s = luaL_checknumber(L, 1);
  lua_Number xn = luaL_checknumber(L, 2);
  lua_Number pr = luaL_checknumber(L, 3);
  lua_Number d;
  check_nbinom(L, 1, s, xn, pr);
  d = exp(xn * log(pr) + s * log(1 - pr) - dlnbet(&s, &xn)) / s;
  lua_pushnumber(L, d);
  return 1;
}

static int stat_pnbinom (lua_State *L) {
  /* stack should contain s, xn, pr */
  lua_Number s = luaL_checknumber(L, 1);
  lua_Number xn = luaL_checknumber(L, 2);
  lua_Number pr = luaL_checknumber(L, 3);
  lua_Number p, q, ompr, bound;
  int which = 1;
  int status;
  check_nbinom(L, 1, s, xn, pr);
  ompr = 1 - pr;
  cdfnbn(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qnbinom (lua_State *L) {
  /* stack should contain p, xn, pr */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number xn = luaL_checknumber(L, 2);
  lua_Number pr = luaL_checknumber(L, 3);
  int si = 0;
  check_nbinom(L, 2, p, xn, pr);
  if (p == 1) {
    lua_pushnumber(L, HUGE_VAL);
    return 1;
  }
  if (p > 0) {
    lua_Number q = 1 - p;
    lua_Number ompr = 1 - pr;
    lua_Number s, bound;
    int which = 2;
    int status;
    cdfnbn(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
    check_status(L, status, bound);
    lua_number2int(si, s);
  }
  lua_pushinteger(L, si);
  return 1;
}


/* {=======   Normal   =======} */

static void check_norm (lua_State *L, int which, lua_Number x,
    lua_Number sd) {
  if (which == 2) luaL_argcheck(L, x >= 0 && x <= 1,  /* p */
      1, "out of range");
  luaL_argcheck(L, sd >= 0, 3, "non-negative value expected");
}

static int stat_dnorm (lua_State *L) {
  /* stack should contain x, and opt. mean and sd */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number mean = luaL_optnumber(L, 2, 0);
  lua_Number sd = luaL_optnumber(L, 3, 1);
  lua_Number d;
  check_norm(L, 1, x, sd);
  d = (x - mean) / sd;
  d = exp(-d*d / 2) / (SQRT2PI * sd);
  lua_pushnumber(L, d);
  return 1;
}

static int stat_pnorm (lua_State *L) {
  /* stack should contain x, and opt. mean and sd */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number mean = luaL_optnumber(L, 2, 0);
  lua_Number sd = luaL_optnumber(L, 3, 1);
  lua_Number p, q, bound;
  int which = 1;
  int status;
  check_norm(L, 1, x, sd);
  q = 1 - p;
  cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qnorm (lua_State *L) {
  /* stack should contain p, and opt. mean and sd */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number mean = luaL_optnumber(L, 2, 0);
  lua_Number sd = luaL_optnumber(L, 3, 1);
  lua_Number x;
  check_norm(L, 2, p, sd);
  if (p == 0 || p == 1) x = (p == 0) ? -HUGE_VAL : HUGE_VAL;
  else {
    lua_Number q = 1 - p;
    lua_Number bound;
    int which = 2;
    int status;
    cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
    check_status(L, status, bound);
  }
  lua_pushnumber(L, x);
  return 1;
}


/* {=======   Poisson   =======} */

static void check_pois (lua_State *L, int which, lua_Number x,
    lua_Number xlam) {
  luaL_argcheck(L, ((which == 1 && x >= 0)  /* x */
      || (which == 2 && (x >= 0 && x <= 1))), /* p */
      1, "out of range");
  luaL_argcheck(L, xlam >= 0, 2, "non-negative value expected");
}

static int stat_dpois (lua_State *L) {
  /* stack should contain s and xlam */
  lua_Number s = luaL_checknumber(L, 1);
  lua_Number xlam = luaL_checknumber(L, 2);
  lua_Number d;
  check_pois(L, 1, s, xlam);
  d = s + 1;
  d = exp(s * log(xlam) - xlam - dlngam(&d));
  lua_pushnumber(L, d);
  return 1;
}

static int stat_ppois (lua_State *L) {
  /* stack should contain s and xlam */
  lua_Number s = luaL_checknumber(L, 1);
  lua_Number xlam = luaL_checknumber(L, 2);
  lua_Number p, q, bound;
  int which = 1;
  int status;
  check_pois(L, 1, s, xlam);
  q = 1 - p;
  cdfpoi(&which, &p, &q, &s, &xlam, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qpois (lua_State *L) {
  /* stack should contain p and xlam */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number xlam = luaL_checknumber(L, 2);
  int si = 0;
  check_pois(L, 2, p, xlam);
  if (p == 1) {
    lua_pushnumber(L, HUGE_VAL);
    return 1;
  }
  if (p > 0) {
    lua_Number q = 1 - p;
    lua_Number s, bound;
    int which = 2;
    int status;
    cdfpoi(&which, &p, &q, &s, &xlam, &status, &bound);
    check_status(L, status, bound);
    lua_number2int(si, s);
  }
  lua_pushinteger(L, si);
  return 1;
}


/* {=======   t-Student   =======} */

static void check_t (lua_State *L, int which, lua_Number x,
    lua_Number df) {
  if (which==2) luaL_argcheck(L, x >= 0 && x <= 1,  /* p */
      1, "out of range");
  luaL_argcheck(L, df >= 0, 3, "non-negative value expected");
}

static int stat_dt (lua_State *L) {
  /* stack should contain x and df */
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number df = luaL_checknumber(L, 2);
  lua_Number t = 0.5;
  lua_Number d;
  check_t(L, 1, x, df);
  d = df / 2;
  d = -dlnbet(&d, &t) - (df + 1) / 2 * log(1 + x * x / df);
  d = exp(d) / sqrt(df);
  lua_pushnumber(L, d);
  return 1;
}

static int stat_pt (lua_State *L) {
  /* stack should contain t and df */
  lua_Number t = luaL_checknumber(L, 1);
  lua_Number df = luaL_checknumber(L, 2);
  lua_Number p, q, bound;
  int which = 1;
  int status;
  check_t(L, 1, t, df);
  q = 1 - p;
  cdft(&which, &p, &q, &t, &df, &status, &bound);
  check_status(L, status, bound);
  lua_pushnumber(L, p);
  return 1;
}

static int stat_qt (lua_State *L) {
  /* stack should contain p and df */
  lua_Number p = luaL_checknumber(L, 1);
  lua_Number df = luaL_checknumber(L, 2);
  lua_Number t;
  check_t(L, 2, p, df);
  if (p == 0 || p == 1) t = (p == 0) ? -HUGE_VAL : HUGE_VAL;
  else {
    lua_Number q = 1 - p;
    lua_Number bound;
    int which = 2;
    int status;
    cdft(&which, &p, &q, &t, &df, &status, &bound);
    check_status(L, status, bound);
  }
  lua_pushnumber(L, t);
  return 1;
}



/* {=====================================================================
 *    Interface
 * ======================================================================} */

static const luaL_Reg factor_lib[] = {
  {"fold", factor_fold},
  {"partition", factor_partition},
  {"design", factor_design},
  {NULL, NULL}
};

static const luaL_Reg stat_lib[] = {
  /* probability dists */
  {"dbeta", stat_dbeta},
  {"pbeta", stat_pbeta},
  {"qbeta", stat_qbeta},
  {"dbinom", stat_dbinom},
  {"pbinom", stat_pbinom},
  {"qbinom", stat_qbinom},
  {"dchisq", stat_dchisq},
  {"pchisq", stat_pchisq},
  {"qchisq", stat_qchisq},
  {"dexp", stat_dexp},
  {"pexp", stat_pexp},
  {"qexp", stat_qexp},
  {"df", stat_df},
  {"pf", stat_pf},
  {"qf", stat_qf},
  {"dgamma", stat_dgamma},
  {"pgamma", stat_pgamma},
  {"qgamma", stat_qgamma},
  {"dhyper", stat_dhyper},
  {"phyper", stat_phyper},
  {"dnbinom", stat_dnbinom},
  {"pnbinom", stat_pnbinom},
  {"qnbinom", stat_qnbinom},
  {"dnorm", stat_dnorm},
  {"pnorm", stat_pnorm},
  {"qnorm", stat_qnorm},
  {"dpois", stat_dpois},
  {"ppois", stat_ppois},
  {"qpois", stat_qpois},
  {"dt", stat_dt},
  {"pt", stat_pt},
  {"qt", stat_qt},
  {NULL, NULL}
};

NUMLUA_API int luaopen_numlua_stat (lua_State *L) {
  luaL_newlib(L, stat_lib);
  /* factors: */
  lua_newtable(L); /* factor mt */
  lua_newtable(L); /* mt */
  lua_pushliteral(L, "k");
  lua_setfield(L, -2, "__mode");
  lua_setmetatable(L, -2); /* factor mt is a weak table */
  /* setup factor mt */
  lua_pushcfunction(L, factor__tostring);
  lua_setfield(L, -2, "__tostring");
  lua_pushcfunction(L, factor__len);
  lua_setfield(L, -2, "__len");
  lua_pushvalue(L, -1); /* factor mt */
  lua_pushcclosure(L, factor__call, 1);
  lua_setfield(L, -2, "__call");
  lua_pushvalue(L, -1); /* factor mt */
  lua_newtable(L);
  nl_register(L, factor_lib, 0);
  lua_pushcclosure(L, factor__index, 2);
  lua_setfield(L, -2, "__index");
  /* setup factor function */
  lua_pushcclosure(L, stat_factor, 1);
  lua_setfield(L, -2, "factor");
  return 1;
}

