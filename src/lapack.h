#ifndef _LAPACK_
#define _LAPACK_

// Diagonalize symmetric matrix
// jobz 'N'= eigenvalues only. 'V' = include eigenvectors
// uplo 'U'= upper triangle storrage, 'L'= lower triangle
// n = matrix order
// a = nxn matrix to diagonalize
// lda = n
// w = the eigenvalues
// work = work array of dimension lwork, work(1) is optimal value of lwork
// lwork dimension of work, if lwork=-1 the optimal value is returned in work(1)
// info 0 = successfull run
extern void ssyev_(
              char *jobz,
              char *uplo,
              int *n,
              float *a,
              int *lda,
              float *w,
              float *work,
              int *lwork,
              int *info
              );

// Diagonalize real nonsymmetric matrix
// [in]	JOBVL	
//           JOBVL is CHARACTER*1
//           = 'N': left eigenvectors of A are not computed;
//           = 'V': left eigenvectors of A are computed.
// [in]	JOBVR	
//           JOBVR is CHARACTER*1
//           = 'N': right eigenvectors of A are not computed;
//           = 'V': right eigenvectors of A are computed.
// [in]	N	
//           N is INTEGER
//           The order of the matrix A. N >= 0.
// [in,out]	A	
//           A is REAL array, dimension (LDA,N)
//           On entry, the N-by-N matrix A.
//           On exit, A has been overwritten.
// [in]	LDA	
//           LDA is INTEGER
//           The leading dimension of the array A.  LDA >= max(1,N).
// [out]	WR	
//           WR is REAL array, dimension (N)
// [out]	WI	
//           WI is REAL array, dimension (N)
//           WR and WI contain the real and imaginary parts,
//           respectively, of the computed eigenvalues.  Complex
//           conjugate pairs of eigenvalues appear consecutively
//           with the eigenvalue having the positive imaginary part
//           first.
// [out]	VL	
//           VL is REAL array, dimension (LDVL,N)
//           If JOBVL = 'V', the left eigenvectors u(j) are stored one
//           after another in the columns of VL, in the same order
//           as their eigenvalues.
//           If JOBVL = 'N', VL is not referenced.
//           If the j-th eigenvalue is real, then u(j) = VL(:,j),
//           the j-th column of VL.
//           If the j-th and (j+1)-st eigenvalues form a complex
//           conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
//           u(j+1) = VL(:,j) - i*VL(:,j+1).
// [in]	LDVL	
//           LDVL is INTEGER
//           The leading dimension of the array VL.  LDVL >= 1; if
//           JOBVL = 'V', LDVL >= N.
// [out]	VR	
//           VR is REAL array, dimension (LDVR,N)
//           If JOBVR = 'V', the right eigenvectors v(j) are stored one
//           after another in the columns of VR, in the same order
//           as their eigenvalues.
//           If JOBVR = 'N', VR is not referenced.
//           If the j-th eigenvalue is real, then v(j) = VR(:,j),
//           the j-th column of VR.
//           If the j-th and (j+1)-st eigenvalues form a complex
//           conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
//           v(j+1) = VR(:,j) - i*VR(:,j+1).
// [in]	LDVR	
//           LDVR is INTEGER
//           The leading dimension of the array VR.  LDVR >= 1; if
//           JOBVR = 'V', LDVR >= N.
// [out]	WORK	
//           WORK is REAL array, dimension (MAX(1,LWORK))
//           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
// [in]	LWORK	
//           LWORK is INTEGER
//           The dimension of the array WORK.  LWORK >= max(1,3*N), and
//           if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
//           performance, LWORK must generally be larger.
//           If LWORK = -1, then a workspace query is assumed; the routine
//           only calculates the optimal size of the WORK array, returns
//           this value as the first entry of the WORK array, and no error
//           message related to LWORK is issued by XERBLA.
// [out]	INFO	
//           INFO is INTEGER
//           = 0:  successful exit
//           < 0:  if INFO = -i, the i-th argument had an illegal value.
//           > 0:  if INFO = i, the QR algorithm failed to compute all the
//                 eigenvalues, and no eigenvectors have been computed;
//                 elements i+1:N of WR and WI contain eigenvalues which
//                 have converged.
extern void sgeev_(
                char 	*JOBVL,
                char 	*JOBVR,
                int 	*N,
                float 	*A,
                int 	*LDA,
                float 	*WR,
                float 	*WI,
                float 	*VL,
                int 	*LDVL,
                float 	*VR,
                int 	*LDVR,
                float 	*WORK,
                int 	*LWORK,
                int 	*INFO     
                );

// Compute LU factorization of a matrix             
// [in]	M	
//           M is INTEGER
//           The number of rows of the matrix A.  M >= 0.
// [in]	N	
//           N is INTEGER
//           The number of columns of the matrix A.  N >= 0.
// [in,out]	A	
//           A is REAL array, dimension (LDA,N)
//           On entry, the M-by-N matrix to be factored.
//           On exit, the factors L and U from the factorization
//           A = P*L*U; the unit diagonal elements of L are not stored.
// [in]	LDA	
//           LDA is INTEGER
//           The leading dimension of the array A.  LDA >= max(1,M).
// [out]	IPIV	
//           IPIV is INTEGER array, dimension (min(M,N))
//           The pivot indices; for 1 <= i <= min(M,N), row i of the
//           matrix was interchanged with row IPIV(i).
// [out]	INFO	
//           INFO is INTEGER
//           = 0:  successful exit
//           < 0:  if INFO = -i, the i-th argument had an illegal value
//           > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
//                 has been completed, but the factor U is exactly
//                 singular, and division by zero will occur if it is used
//                 to solve a system of equations.   
extern void sgetrf_(
                int* 	M,
                int* 	N,
                float* 	A,
                int* 	LDA,
                int* 	IPIV,
                int* 	INFO 
                );

// Compute the inverse of a matrix using LU factorization
// [in]	N	
//           N is INTEGER
//           The order of the matrix A.  N >= 0.
// [in,out]	A	
//           A is REAL array, dimension (LDA,N)
//           On entry, the factors L and U from the factorization
//           A = P*L*U as computed by SGETRF.
//           On exit, if INFO = 0, the inverse of the original matrix A.
// [in]	LDA	
//           LDA is INTEGER
//           The leading dimension of the array A.  LDA >= max(1,N).
// [in]	IPIV	
//           IPIV is INTEGER array, dimension (N)
//           The pivot indices from SGETRF; for 1<=i<=N, row i of the
//           matrix was interchanged with row IPIV(i).
// [out]	WORK	
//           WORK is REAL array, dimension (MAX(1,LWORK))
//           On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
// [in]	LWORK	
//           LWORK is INTEGER
//           The dimension of the array WORK.  LWORK >= max(1,N).
//           For optimal performance LWORK >= N*NB, where NB is
//           the optimal blocksize returned by ILAENV.
//           If LWORK = -1, then a workspace query is assumed; the routine
//           only calculates the optimal size of the WORK array, returns
//           this value as the first entry of the WORK array, and no error
//           message related to LWORK is issued by XERBLA.
// [out]	INFO	
//           INFO is INTEGER
//           = 0:  successful exit
//           < 0:  if INFO = -i, the i-th argument had an illegal value
//           > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
//                 singular and its inverse could not be computed.
extern void sgetri_(
                int* 	N,
                float* 	A,
                int* 	LDA,
                int* 	IPIV,
                float* 	WORK,
                int* 	LWORK,
                int* 	INFO 
                );

#endif // LAPACK