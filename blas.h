/* ============================= Level 1 ============================= */

double ddot_(int* N, double* X, int* incX, double* Y,
    int* incY);

void zdotu_(double complex* retval, int* N, double complex* X, int*
    incX, double complex* Y, int* incY);

void zdotc_(double complex* retval, int* N, double complex* X, int*
    incX, double complex* Y, int* incY);

double dnrm2_(int* N, double* X, int* incX);

double dasum_(int* N, double* X, int* incX);

double dznrm2_(int* N, double complex* X, int* incX);

double dzasum_(int* N, double complex* X, int* incX);

int idamax_(int* N, double* X, int* incX);

int izamax_(int* N, double complex* X, int* incX);

int dswap_(int* N, double* X, int* incX, double* Y, int*
    incY);

int dcopy_(int* N, double* X, int* incX, double* Y, int*
    incY);

int daxpy_(int* N, double* alpha, double* X, int* incX,
    double* Y, int* incY);

int zswap_(int* N, double complex* X, int* incX, double complex* Y,
    int* incY);

int zcopy_(int* N, double complex* X, int* incX, double complex* Y,
    int* incY);

int zaxpy_(int* N, double complex* alpha, double complex* X, int* incX,
    double complex* Y, int* incY);

int drotg_(double* a, double* b, double* c, double* s);

int drot_(int* N, double* X, int* incX, double* Y, int*
    incY, double* c, double* s);

int dscal_(int* N, double* alpha, double* X, int* incX);

int zscal_(int* N, double complex* alpha, double complex* X, int* incX);

int zdscal_(int* N, double* alpha, double complex* X, int* incX);

/* ============================= Level 2 ============================= */

int dgemv_(char* trans, int* M, int* N, double* alpha, double*
    A, int* lda, double* X, int* incX, double* beta,
    double* Y, int* incY, int ltrans);

int dtrmv_(char* uplo, char *trans, char* diag, int *N,  double *A,
    int *lda, double *X, int *incX, int luplo, int ltrans,
    int ldiag);

int dtrsv_(char* uplo, char* trans, char* diag, int* N, double* A,
    int* lda, double* X, int* incX, int luplo, int ltrans,
    int ldiag);

int zgemv_(char* trans, int* M, int* N, double complex* alpha,
    double complex* A, int* lda, double complex* X, int* incX,
    double complex* beta, double complex* Y, int* incY, int ltrans);

int ztrmv_(char* uplo, char *trans, char* diag, int *N,  double complex *A,
    int *lda, double complex *X, int *incX, int luplo, int ltrans,
    int ldiag);

int ztrsv_(char* uplo, char* trans, char* diag, int* N, double complex* A,
    int* lda, double complex* X, int* incX, int luplo, int ltrans,
    int ldiag);

int dsymv_(char* uplo, int* N, double* alpha, double* A, int*
    lda, double* X, int* incX, double* beta, double* Y,
    int* incY, int luplo);

int dger_(int* M, int* N, double* alpha, double* X, int*
    incX, double* Y, int* incY, double* A, int* lda);

int dsyr_(char* uplo, int* N, double* alpha, double* X, int*
    incX, double* A, int* lda, int luplo);

int dsyr2_(char* uplo, int* N, double* alpha, double* X, int*
    incX, double* Y, int* incY, double* A, int* lda,
    int luplo);

int zhemv_(char* uplo, int* N, double complex* alpha, double complex* A,
    int* lda, double complex* X, int* incX, double complex* beta,
    double complex* Y, int* incY, int luplo);

int zgeru_(int* M, int* N, double complex* alpha, double complex* X,
    int* incX, double complex* Y, int* incY, double complex* A,
    int* lda);

int zgerc_(int* M, int* N, double complex* alpha, double complex* X,
    int* incX, double complex* Y, int* incY, double complex* A,
    int* lda);

int zher_(char* uplo, int* N, double* alpha, double complex* X,
    int* incX, double complex* A, int* lda, int luplo);

int zher2_(char* uplo, int* N, double complex* alpha, double complex* X,
    int* incX, double complex* Y, int* incY, double complex* A,
    int* lda, int luplo);

/* ============================= Level 3 ============================= */

int dgemm_(char* transA, char* transB, int* M, int* N, int* K,
    double* alpha, double* A, int* lda, double* B,
    int* ldb, double* beta, double* C, int* ldc,
    int ltransa, int ltransb);

int dsymm_(char* side, char* uplo, int* M, int* N, double* alpha,
    double* A, int* lda, double* B, int* ldb, double*
    beta, double* C, int* ldc, int lside, int luplo);

int dsyrk_(char* uplo, char* trans, int* N, int* K, double* alpha,
    double* A, int* lda, double* beta, double* C, int*
    ldc, int luplo, int ltrans);

int dsyr2k_(char* uplo, char* trans, int* N, int* K, double*
    alpha, double* A, int* lda, double* B, int* ldb,
    double* beta, double* C, int* ldc, int luplo,
    int ltrans);

int dtrmm_(char* side, char* uplo, char* trans, char* diag, int* M,
    int* N, double* alpha, double* A, int* lda,
    double* B, int* ldb, int lside, int luplo,
    int ltrans, int ldiag);

int dtrsm_(char* side, char* uplo, char* trans, char* diag, int* M,
    int* N, double* alpha, double* A, int* lda,
    double* B, int* ldb, int lside, int luplo,
    int ltrans, int ldiag);

int zgemm_(char* transA, char* transB, int* M, int* N, int* K,
    double complex* alpha, double complex* A, int* lda, double complex*
    B, int* ldb, double complex* beta, double complex* C, int* ldc,
    int ltransa, int ltransb);

int zsyrk_(char* uplo, char* trans, int* N, int* K, double complex*
    alpha, double complex* A, int* lda, double complex* beta,
    double complex* C, int* ldc, int luplo, int ltrans);

int zsyr2k_(char* uplo, char* trans, int* N, int* K, double complex*
    alpha, double complex* A, int* lda, double complex* B, int* ldb,
    double complex* beta, double complex* C, int* ldc, int luplo,
    int ltrans);

int ztrmm_(char* side, char* uplo, char* trans, char* diag, int* M,
    int* N, double complex* alpha, double complex* A, int* lda,
    double complex* B, int* ldb, int lside, int luplo,
    int ltrans, int ldiag);

int ztrsm_(char* side, char* uplo, char* trans, char* diag, int* M,
    int* N, double complex* alpha, double complex* A, int* lda,
    double complex* B, int* ldb, int lside, int luplo,
    int ltrans, int ldiag);

int zhemm_(char* side, char* uplo, int* M, int* N, double complex*
    alpha, double complex* A, int* lda, double complex* B, int* ldb,
    double complex* beta, double complex* C, int* ldc, int lside,
    int luplo);

int zherk_(char* uplo, char* trans, int* N, int* K, double* alpha,
    double complex* A, int* lda, double* beta, double complex* C,
    int* ldc, int luplo, int ltrans);

int zher2k_(char* uplo, char* trans, int* N, int* K, double complex*
    alpha, double complex* A, int* lda, double complex* B, int* ldb,
    double* beta, double complex* C, int* ldc, int luplo,
    int ltrans);

