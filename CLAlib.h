#define F77NAME(x) x##_

extern "C"
{
	void F77NAME(dgemv)(const char &trans, const int &m, const int &n,
			const double &alpha, const double *A, const int &lda, 
			const double *x, const int &incx, const double &beta,
			const double *y, const int &incy);
	void F77NAME(daxpy) (const int &n, const double &alpha,
			const double *x, const int &incx,
		       	const double *y, const int &incy);
}
