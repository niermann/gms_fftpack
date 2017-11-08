#ifndef FFTPACK_FFTPACK_INC
#define FFTPACK_FFTPACK_INC

// restrict keyword is only supported since MSVC 2005 */
#if defined(_MSC_VER) && (_MSC_VER < 1400)
#       define RESTRICT
#else
#       define RESTRICT __restrict
#endif

/**
 * subroutine cfftf computes the forward complex discrete fourier
 * transform (the fourier analysis). equivalently , cfftf computes
 * the fourier coefficients of a complex periodic sequence.
 * the transform is defined below at output parameter c.
 * 
 * the transform is not normalized. to obtain a normalized transform
 * the output must be divided by n. otherwise a call of cfftf
 * followed by a call of cfftb will multiply the sequence by n.
 * 
 * the array wsave which is used by subroutine cfftf must be
 * initialized by calling subroutine cffti(n,wsave).
 * 
 * input parameters
 * n       the length of the complex sequence c. the method is
 *         more efficient when n is the product of small primes. n
 * c       a complex array of length n which contains the sequence
 * ch      a real work array which must be dimensioned at least 2n and 
 *         will be used for temporary storage during the call
 * wsave   a real array which must be dimensioned at least 2n
 * facsave an int array which must be dimensioned at least 15 elements
 *         The wsave and facsave arrays must be initialized by calling 
 *         subroutine cffti(n,wsave,facsave) and a different wsave/facsave array 
 *         must be used for each different value of n. this initialization 
 *         does not have to be repeated so long as n remains unchanged thus 
 *         subsequent transforms can be obtained faster than the first.
 *         The same wsave/facsave array can be used by cfftf and cfftb.
 * 
 * output parameters
 * c       c(j)=sum(0 <= k < n) of
 *                c(k)*exp(-i*j*k*2*pi/n)
 *         where i=sqrt(-1)
 **/
extern void cfftf(int n, float* RESTRICT c, float* RESTRICT ch, const float* RESTRICT wsave, const int* RESTRICT facsave);
extern void dcfftf(int n, double* RESTRICT c, double* RESTRICT ch, const double* RESTRICT wsave, const int* RESTRICT facsave);

/**
 * subroutine cfftb computes the backward complex discrete fourier
 * transform (the fourier synthesis). equivalently , cfftb computes
 * a complex periodic sequence from its fourier coefficients.
 * the transform is defined below at output parameter c.
 * 
 * a call of cfftf followed by a call of cfftb will multiply the
 * sequence by n.
 * 
 * the array wsave which is used by subroutine cfftb must be
 * initialized by calling subroutine cffti(n,wsave).
 *
 * input parameters
 * n       the length of the complex sequence c. the method is
 *         more efficient when n is the product of small primes.
 *         a complex array of length n which contains the sequence
 * ch      a real work array which must be dimensioned at least 2n and 
 *         will be used for temporary storage during the call
 * wsave   a real array which must be dimensioned at least 2n
 * facsave an int array which must be dimensioned at least 15 elements
 *         The wsave and facsave arrays must be initialized by calling 
 *         subroutine cffti(n,wsave,facsave) and a different wsave/facsave array 
 *         must be used for each different value of n. this initialization 
 *         does not have to be repeated so long as n remains unchanged thus 
 *         subsequent transforms can be obtained faster than the first.
 *         The same wsave/facsave array can be used by cfftf and cfftb.
 *
 * output parameters
 * c       c(j)=sum(0 <= k < n) of
 *                c(k)*exp(+i*j*k*2*pi/n)
 *         where i=sqrt(-1)
 **/
extern void cfftb(int n, float* RESTRICT c, float* RESTRICT ch, const float* RESTRICT wsave, const int* RESTRICT facsave);
extern void dcfftb(int n, double* RESTRICT c, double* RESTRICT ch, const double* RESTRICT wsave, const int* RESTRICT facsave);

/**
 * subroutine cffti initializes the array wsave which is used in
 * both cfftf and cfftb. the prime factorization of n together with
 * a tabulation of the trigonometric functions are computed and
 * stored in wsave.
 * 
 * input parameter
 * n       the length of the sequence to be transformed
 * 
 * output parameter
 * wsave   a array which must be dimensioned at least 2*n
 * facsave a array which must be dimensioned at least 15
 *         the same rrays can be used for both cfftf and cfftb
 *         as long as n remains unchanged. different wsave/facsave arrays
 *         are required for different values of n. the contents of
 *         wsave must not be changed between calls of cfftf or cfftb.
 **/
extern void cffti(int n, float* RESTRICT wsave, int* RESTRICT facsave);
extern void dcffti(int n, double* RESTRICT wsave, int* RESTRICT facsave);

/**
 * subroutine rfftf computes the fourier coefficients of a real
 * perodic sequence (fourier analysis). the transform is defined
 * below at output parameter r.
 * 
 * input parameters
 * n       the length of the array r to be transformed.  the method
 *         is most efficient when n is a product of small primes.
 *         n may change so long as different work arrays are provided
 * r       a real array of length n which contains the sequence
 *         to be transformed
 * ch      a real work array which must be dimensioned at least n and 
 *         will be used for temporary storage during the call
 * wsave   a real array which must be dimensioned at least n
 * facsave an int array which must be dimensioned at least 15 elements
 *         The wsave and facsave arrays must be initialized by calling 
 *         subroutine rffti(n,wsave,facsave) and a different wsave/facsave array 
 *         must be used for each different value of n. this initialization 
 *         does not have to be repeated so long as n remains unchanged thus 
 *         subsequent transforms can be obtained faster than the first.
 *         The same wsave/facsave array can be used by rfftf and rfftb.
 *
 * output parameters
 * r       r(0) = sum(0 <= j < n) of r(j)
 *         if n is even set L=n/2, if n is odd set L=(n+1)/2
 *         then for 1 <= k < L
 *           r(2*k-1) = sum(0 <= j < n) of
 *                r(j)*cos(k*j*2*pi/n)
 *           r(2*k) = sum(0 <= j < n) of
 *               -r(i)*sin(k*j*2*pi/n)
 *         if n is even
 *           r(n) = sum(0 <= j < n) of 
 *                (-1)**j*r(j)
 *
 * note
 *         this transform is unnormalized since a call of rfftf
 *         followed by a call of rfftb will multiply the input
 *         sequence by n.
 **/
extern void rfftf(int n, float* r, float* ch, const float* RESTRICT wasave, const int* RESTRICT facsave);
extern void drfftf(int n, double* r, double* ch, const double* RESTRICT wasave, const int* RESTRICT facsave);

/**
 * subroutine rfftb computes the real perodic sequence from its
 * fourier coefficients (fourier synthesis). the transform is defined
 * below at output parameter r.
 * 
 * input parameters
 * n       the length of the array r to be transformed.  the method
 *         is most efficient when n is a product of small primes.
 *         n may change so long as different work arrays are provided
 * c       a real array of length n which contains the sequence
 *         to be transformed
 * ch      a real work array which must be dimensioned at least n and 
 *         will be used for temporary storage during the call
 * wsave   a real array which must be dimensioned at least n
 * facsave an int array which must be dimensioned at least 15 elements
 *         The wsave and facsave arrays must be initialized by calling 
 *         subroutine rffti(n,wsave,facsave) and a different wsave/facsave array 
 *         must be used for each different value of n. this initialization 
 *         does not have to be repeated so long as n remains unchanged thus 
 *         subsequent transforms can be obtained faster than the first.
 *         The same wsave/facsave array can be used by rfftf and rfftb.
 *
 * output parameters
 * r       for n even
 *             r(j) = r(0) + (-1)**j*r(n-1)
 *                  + sum(1 <= k < n/2) of
 *                         2.*r(2*k)*cos(k*j*2*pi/n)
 *                        -2.*r(2*k+1)*sin(k*j*2*pi/n)
 *         for n odd
 *             r(j) = r(0) 
 *                  + sum(1 <= k < (n+1)/2) of
 *                        2.*r(2*k)*cos(k*j*2*pi/n)
 *                       -2.*r(2*k+1)*sin(k*j*2*pi/n)
 *
 * note
 *         this transform is unnormalized since a call of rfftf
 *         followed by a call of rfftb will multiply the input
 *         sequence by n.
 **/
extern void rfftb(int n, float* r, float* ch, const float* RESTRICT wasave, const int* RESTRICT facsave);
extern void drfftb(int n, double* r, double* ch, const double* RESTRICT wasave, const int* RESTRICT facsave);

/**
 * subroutine rffti initializes the array wsave which is used in
 * both rfftf and rfftb. the prime factorization of n together with
 * a tabulation of the trigonometric functions are computed and
 * stored in wsave.
 * 
 * input parameter
 * n       the length of the sequence to be transformed.
 * 
 * output parameter
 * wsave   a real array which must be dimensioned at least n
 * facsave an int array which must be dimensioned at least 15 elements
 *         the same arrays can be used for both rfftf and rfftb
 *         as long as n remains unchanged. different wsave/facsave arrays
 *         are required for different values of n. the contents of
 *         wsave must not be changed between calls of rfftf or rfftb.
 **/
extern void rffti(int n, float* RESTRICT wasave, int* RESTRICT facsave);
extern void drffti(int n, double* RESTRICT wasave, int* RESTRICT facsave);

#endif // FFTPACK_FFTPACK_INC
