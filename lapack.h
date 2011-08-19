/* =========================== Internals ============================= */

/* DLAMCH - determine double precision machine parameters */
double dlamch_(char *cmach, int lmach);

/* xLANGE - return the value of the one norm, or the Frobenius norm, or the
 * infinity norm, or the element of largest absolute value of a real matrix A
 * */
double dlange_(char *norm, int *m, int *n, double *a,
    int *lda, double *work, int lnorm);
double zlange_(char *norm, int *m, int *n, double complex *a,
    int *lda, double *work, int lnorm);

/* xLASET - initialize an m-by-n matrix A to BETA on the diagonal and
 * ALPHA on the offdiagonals */
int dlaset_(char *uplo, int *m, int *n, double *alpha,
    double *beta, double *a, int *lda, int luplo);
int zlaset_(char *uplo, int *m, int *n, double complex *alpha,
    double complex *beta, double complex *a, int *lda, int luplo);

/* =========================== General ============================= */

/* xGEBAK - form the right or left eigenvectors of a real [complex] general
 * matrix by backward transformation on the computed eigenvectors of the
 * balanced matrix output by xGEBAL */
int dgebak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, double *scale, int *m, double *v,
    int *ldv, int *info, int ljob, int lside);

int zgebak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, double *scale, int *m, double complex *v, 
    int *ldv, int *info, int ljob, int lside);

/* xGEBAL - balance a general real [complex] matrix A */
int dgebal_(char *job, int *n, double *a, int *lda,
    int *ilo, int *ihi, double *scale, int *info,
    int ljob);

int zgebal_(char *job, int *n, double complex *a, int *lda,
    int *ilo, int *ihi, double *scale, int *info,
    int ljob);

/* xGECON - estimate the reciprocal of the condition number of a general real
 * [complex] matrix A, in either the 1-norm or the infinity-norm, using the LU
 * factorization computed by xGETRF */
int dgecon_(char *norm, int *n, double *a, int *lda,
    double *anorm, double *rcond, double *work,
    int *iwork, int *info, int lnorm);

int zgecon_(char *norm, int *n, double complex *a, int *lda,
    double *anorm, double *rcond, double complex *work,
    double *rwork, int *info, int lnorm); 

/* xGEEV - compute for an N-by-N real [complex] nonsymmetric matrix A, the
 * eigenvalues and, optionally, the left and/or right eigenvectors */
int dgeev_(char *jobvl, char *jobvr, int *n, double *a,
    int *lda, double *wr, double *wi, double *vl, 
    int *ldvl, double *vr, int *ldvr, double *work, 
    int *lwork, int *info, int ljobvl, int ljobvr);

int zgeev_(char *jobvl, char *jobvr, int *n, double complex *a,
    int *lda, double complex *w, double complex *vl, int *ldvl,
    double complex *vr, int *ldvr, double complex *work, 
    int *lwork, double *rwork, int *info, int ljobvl,
    int ljobvr);
 
/* xGEHRD - reduce a real [complex] general matrix A to upper Hessenberg form
 * H by an orthogonal similarity transformation */
int dgehrd_(int *n, int *ilo, int *ihi, double *a,
    int *lda, double *tau, double *work, 
    int *lwork, int *info);

int zgehrd_(int *n, int *ilo, int *ihi, double complex *a,
    int *lda, double complex *tau, double complex *work,
    int *lwork, int *info);
 
/* xGELQF - compute an LQ factorization of a real [complex] M-by-N matrix A */
int dgelqf_(int *m, int *n, double *a, int *lda,
    double *tau, double *work, int *lwork, int *info);

int zgelqf_(int *m, int *n, double complex *a, int *lda,
    double complex *tau, double complex *work, int *lwork,
    int *info);
 
/* xGELS - solve overdetermined or underdetermined real [complex] linear
 * systems involving an M-by-N matrix A, or its transpose, using a QR or LQ
 * factorization of A */
int dgels_(char *trans, int *m, int *n, int *nrhs,
    double *a, int *lda, double *b, int *ldb, 
    double *work, int *lwork, int *info, int ltrans);

int zgels_(char *trans, int *m, int *n, int *nrhs,
    double complex *a, int *lda, double complex *b, int *ldb, 
    double complex *work, int *lwork, int *info, int ltrans);

/* xGELSD - compute the minimum-norm solution to a real [complex] linear
 * least squares problem */
int dgelsd_(int *m, int *n, int *nrhs, double *a,
    int *lda, double *b, int *ldb, double *s,
    double *rcond, int *rank, double *work, int *lwork,
    int *iwork, int *info);

int zgelsd_(int *m, int *n, int *nrhs, double complex *a,
    int *lda, double complex *b, int *ldb, double *s,
    double *rcond, int *rank, double complex *work, int *lwork,
    double *rwork, int *iwork, int *info);

/* xGELSY  - compute the minimum-norm solution to a real [complex] linear
 * least squares problem using complete orthogonal decomposition */
int dgelsy_(int *m, int *n, int *nrhs, double *a,
    int *lda, double *b, int *ldb, int *jpvt,
    double *rcond, int *rank, double *work, int *lwork,
    int *info);

int zgelsy_(int *m, int *n, int *nrhs, double complex *a,
    int *lda, double complex *b, int *ldb, int *jpvt,
    double *rcond, int *rank, double complex *work, int *lwork,
    double *rwork, int *info);

/* xGELSS  - compute the minimum norm solution to a real [complex] linear
 * least squares problem using SVD decomposition */
int dgelss_(int *m, int *n, int *nrhs, double *a,
    int *lda, double *b, int *ldb, double *s,
    double *rcond, int *rank, double *work, int *lwork,
    int *info);

int zgelss_(int *m, int *n, int *nrhs, double complex *a,
    int *lda, double complex *b, int *ldb, double *s,
    double *rcond, int *rank, double complex *work, int *lwork,
    double *rwork, int *info);

/* xGEQP3 - compute a QR factorization with column pivoting of a matrix A */
int dgeqp3_(int *m, int *n, double *a, int *lda,
    int *jpvt, double *tau, double *work, int *lwork,
    int *info);

int zgeqp3_(int *m, int *n, double complex *a, int *lda,
    int *jpvt, double complex *tau, double complex *work, 
    int *lwork, double *rwork, int *info);

/* xGEQRF - compute a QR factorization of a real [complex] M-by-N matrix A */
int dgeqrf_(int *m, int *n, double *a, int *lda,
    double *tau, double *work, int *lwork, int *info);

int zgeqrf_(int *m, int *n, double complex *a, int *lda,
    double complex *tau, double complex *work, int *lwork,
    int *info);
 
/* xGESV - compute the solution to a real [complex] system of linear equations
 * A * X = B, where A is an N-by-N matrix and X and B are N-by-NRHS matrices */
int dgesv_(int *n, int *nrhs, double *a, int *lda,
    int *ipiv, double *b, int *ldb, int *info);

int zgesv_(int *n, int *nrhs, double complex *a, int *lda,
    int *ipiv, double complex *b, int *ldb, int *info);
 
/* xGESVD - compute the singular value decomposition (SVD) of a real [complex]
 * M-by-N matrix A, optionally computing the left and/or right singular
 * vectors */
int dgesvd_(char *jobu, char *jobvt, int *m, int *n, 
    double *a, int *lda, double *s, double *u,
    int *ldu, double *vt, int *ldvt, double *work,
    int *lwork, int *info, int ljobu, int ljobvt);

int zgesvd_(char *jobu, char *jobvt, int *m, int *n, 
    double complex *a, int *lda, double *s, double complex *u,
    int *ldu, double complex *vt, int *ldvt, double complex *work,
    int *lwork, double *rwork, int *info, int ljobu,
    int ljobvt);

/* xGETRF - compute an LU factorization of a general M-by-N matrix A using
 * partial pivoting with row interchanges */
int dgetrf_(int *m, int *n, double *a, int *lda,
    int *ipiv, int *info);

int zgetrf_(int *m, int *n, double complex *a, int *lda,
    int *ipiv, int *info);

/* xGETRI - compute the inverse of a matrix using the LU factorization
 * computed by xGETRF */
int dgetri_(int *n, double *a, int *lda, int *ipiv,
    double *work, int *lwork, int *info);

int zgetri_(int *n, double complex *a, int *lda, int *ipiv,
    double complex *work, int *lwork, int *info);

/* xGETRS - solve a system of linear equations A * X = B or A' * X = B with a
 * general N-by-N matrix A using the LU factorization computed by xGETRF */
int dgetrs_(char *trans, int *n, int *nrhs, double *a,
    int *lda, int *ipiv, double *b, int *ldb,
    int *info, int ltrans);

int zgetrs_(char *trans, int *n, int *nrhs, double complex *a,
    int *lda, int *ipiv, double complex *b, int *ldb,
    int *info, int ltrans);
 
/* DORGHR - generate a real orthogonal matrix Q which is defined as the
 * product of IHI-ILO elementary reflectors of order N, as returned by DGEHRD
 * */
int dorghr_(int *n, int *ilo, int *ihi, double *a,
    int *lda, double *tau, double *work, int *lwork,
    int *info);

/* DORGLQ - generate an M-by-N real matrix Q with orthonormal rows, which is
 * defined as the first M rows of a product of K elementary reflectors of
 * order N as returned by DGELQF */
int dorglq_(int *m, int *n, int *k, double *a,
    int *lda, double *tau, double *work, int *lwork, 
    int *info);

/* xOR[UN]GQR - generate an M-by-N real [complex] matrix Q with orthonormal
 * columns, which is defined as the first N columns of a product of K
 * elementary reflectors of order M as returned by xGEQRF */
int dorgqr_(int *m, int *n, int *k, double *a,
    int *lda, double *tau, double *work, int *lwork, 
    int *info);

int zungqr_(int *m, int *n, int *k, double complex *a,
    int *lda, double complex *tau, double complex *work, int *lwork, 
    int *info);

/* DORMHR - overwrite the general real M-by-N matrix C with SIDE = 'L' SIDE =
 * 'R' TRANS = 'N': Q * C C * Q TRANS = 'T': Q**T * C    C * Q**T where Q is a
 * real orthogonal matrix of order nq, with nq = m if SIDE = 'L' and nq = n if
 * SIDE = 'R'. Q is defined as the product of IHI-ILO elementary reflectors,
 * as returned by DGEHRD */
int dormhr_(char *side, char *trans, int *m, int *n, 
  int *ilo, int *ihi, double *a, int *lda, double *
  tau, double *c, int *ldc, double *work, int *lwork, 
  int *info, int lside, int ltrans);

/* DORMLQ - overwrite the general real M-by-N matrix C with SIDE = 'L' SIDE =
 * 'R' TRANS = 'N': Q * C C * Q TRANS = 'T': Q**T * C   C * Q**T where Q is a
 * real orthogonal matrix defined as the product of k elementary reflectors,
 * as returned by DGELQF */
int dormlq_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau,
    double *c, int *ldc, double *work, int *lwork,
    int *info, int lside, int ltrans);

/* xOR[UN]MQR - overwrite the general real M-by-N matrix C with SIDE = 'L' SIDE =
 * 'R' TRANS = 'N': Q * C C * Q TRANS = 'T': Q**T * C   C * Q**T where Q is a
 * real orthogonal matrix defined as the product of k elementary reflectors,
 * as returned by xGEQRF */
int dormqr_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau,
    double *c, int *ldc, double *work, int *lwork,
    int *info, int lside, int ltrans);

int zunmqr_(char *side, char *trans, int *m, int *n, 
    int *k, double complex *a, int *lda, double complex *tau,
    double complex *c, int *ldc, double complex *work, int *lwork,
    int *info, int lside, int ltrans);


/* =========================== Posdef Symm ============================= */

/* xPOCON - estimate the reciprocal of the condition number (in the 1-norm) of
 * a real [complex] symmetric [Hermitian] positive definite matrix using the
 * Cholesky factorization A = U**H*U or A = L*L**H computed by xPOTRF */
int dpocon_(char *uplo, int *n, double *a, int *lda,
    double *anorm, double *rcond, double *work,
    int *iwork, int *info, int luplo);

int zpocon_(char *uplo, int *n, double complex *a, int *lda,
    double *anorm, double *rcond, double complex *work,
    double *rwork, int *info, int luplo);
 
/* xPOTRF - compute the Cholesky factorization of a real [complex] symmetric
 * [Hermitian] positive definite matrix A */
int dpotrf_(char *uplo, int *n, double *a, int *lda,
    int *info, int luplo);
 
int zpotrf_(char *uplo, int *n, double complex *a, int *lda,
    int *info, int luplo);

/* xPOTRI - compute  the  inverse of a real [complex] symmetric [Hermitian]
 * positive definite matrix A using the Cholesky factorization A = U**H*U or
 * A  = L*L**H computed by xPOTRF */
int dpotri_(char *uplo, int *n, double *a, int *lda,
    int *info, int luplo);

int zpotri_(char *uplo, int *n, double complex *a, int *lda,
    int *info, int luplo);
 
/* xPOTRS - solve a system of linear equations A*X = B with a symmetric
 * [Hermitian] positive definite matrix A using the Cholesky factorization
 * A = U**H*U or A = L*L**H computed by xPOTRF */
int dpotrs_(char *uplo, int *n, int *nrhs, double *a,
    int *lda, double *b, int *ldb, int *info, int luplo);

int zpotrs_(char *uplo, int *n, int *nrhs, double complex *a,
    int *lda, double complex *b, int *ldb, int *info, int luplo);


/* =========================== Symmetric ============================= */

/* xSY[HE]CON - estimate the reciprocal of the condition number (in the 1-norm)
 * of a real symmetric [SY] or complex Hermitian [HE] matrix A using the
 * factorization A = U*D*U**H or A = L*D*L**H computed by xSY[HE]TRF */
int dsycon_(char *uplo, int *n, double *a, int *lda,
    int *ipiv, double *anorm, double *rcond, double *work,
    int *iwork, int *info, int luplo);

int zhecon_(char *uplo, int *n, double complex *a, int *lda,
    int *ipiv, double *anorm, double *rcond, 
    double complex *work, int *info, int luplo);
 
/* xSY[HE]EV - compute all eigenvalues and, optionally, eigenvectors of a real
 * symmetric [SY] or complex Hermitian [HE] matrix A */
int dsyev_(char *jobz, char *uplo, int *n, double *a,
    int *lda, double *w, double *work, int *lwork,
    int *info, int ljobz, int luplo);

int zheev_(char *jobz, char *uplo, int *n, double complex *a,
    int *lda, double *w, double complex *work, int *lwork, 
    double *rwork, int *info, int ljobz, int luplo);

/* xSY[HE]SV - compute the solution to a real [complex] system of linear
 * equations A * X = B, where A is an N-by-N symmetric [Hermitian] matrix and
 * X and B are N-by-NRHS matrices */
int dsysv_(char *uplo, int *n, int *nrhs, double *a,
    int *lda, int *ipiv, double *b, int *ldb, 
    double *work, int *lwork, int *info, int luplo);

int zhesv_(char *uplo, int *n, int *nrhs, double complex *a,
    int *lda, int *ipiv, double complex *b, int *ldb,
    double complex *work, int *lwork, int *info, int luplo);
 
/* xSY[HE]TRF - compute the factorization of a real symmetric [SY] or complex
 * Hermitian [HE] matrix A using the Bunch-Kaufman diagonal pivoting method */
int dsytrf_(char *uplo, int *n, double *a, int *lda,
    int *ipiv, double *work, int *lwork, int *info,
    int luplo);

int zhetrf_(char *uplo, int *n, double complex *a, int *lda,
    int *ipiv, double complex *work, int *lwork, int *info,
    int luplo);
 
/* xSY[HE]TRI - compute the inverse of a real symmetric [SY] or complex
 * Hermitian [HE] indefinite matrix A using the factorization A = U*D*U**H
 * or A = L*D*L**H computed by xSY[HE]TRF */
int dsytri_(char *uplo, int *n, double *a, int *lda,
    int *ipiv, double *work, int *info, int luplo);

int zhetri_(char *uplo, int *n, double complex *a, int *lda,
    int *ipiv, double complex *work, int *info, int luplo);
 
/* xSY[HE]TRS - solve a system of linear equations A*X = B with a real
 * symmetric [SY] or complex Hermitian [HE] matrix A using the factorization
 * A = U*D*U**H or A = L*D*L**H computed by xSY[HE]TRF */
int dsytrs_(char *uplo, int *n, int *nrhs, double *a,
    int *lda, int *ipiv, double *b, int *ldb,
    int *info, int luplo);

int zhetrs_(char *uplo, int *n, int *nrhs, double complex *a,
    int *lda, int *ipiv, double complex *b, int *ldb,
    int *info, int luplo);
 

/* =========================== Triangular ============================= */

/* xTRCON - estimate the reciprocal of the condition number of a triangular
 * matrix A, in either the 1-norm or the infinity-norm */
int dtrcon_(char *norm, char *uplo, char *diag, int *n, 
    double *a, int *lda, double *rcond, double *work, 
    int *iwork, int *info, int lnorm, int luplo,
    int ldiag);

int ztrcon_(char *norm, char *uplo, char *diag, int *n, 
    double complex *a, int *lda, double *rcond,
    double complex *work, double *rwork, int *info, int lnorm,
    int luplo, int ldiag);
 
/* xTRTRI - compute the inverse of a real [complex] upper or lower triangular
 * matrix A */
int dtrtri_(char *uplo, char *diag, int *n, double *a,
    int *lda, int *info, int luplo, int ldiag);

int ztrtri_(char *uplo, char *diag, int *n, double complex *a,
    int *lda, int *info, int luplo, int ldiag);

/* xTRTRS - solve a triangular system of the form A * X = B or A**H * X = B,
 * where A is a triangular matrix of order N, and B is an N-by-NRHS matrix. A
 * check is made to verify that A is nonsingular */
int dtrtrs_(char *uplo, char *trans, char *diag, int *n,
    int *nrhs, double *a, int *lda, double *b,
    int *ldb, int *info, int luplo, int ltrans,
    int ldiag);

int ztrtrs_(char *uplo, char *trans, char *diag, int *n,
    int *nrhs, double complex *a, int *lda, double complex *b, 
    int *ldb, int *info, int luplo, int ltrans,
    int ldiag);
 
