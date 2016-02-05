////////////////////////////////////////////////////////////////
// Filename: fft.c
//
// Synopsis: Example FFT routine. Can run on nearly any CPU that
//   supports a C compiler. This code calculates the FFT
//   (decimation-in-time) of an N point cfloat data sequence.

////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>

#define PI 3.1415926535897932


typedef struct complex_float {
    float re;  // The Real part of a complex number
    float im;  // The Imaginary part of a complex number
} cfloat;

////////////////////////////////////////////


void fft_c(int n, cfloat *x, cfloat *W)
///////////////////////////////////////////////////////////////////////
// Purpose:   Calculate the radix-2 decimation-in-frequency FFT.
//
// Input:     n: length of FFT, x: input array of cfloat numbers,
//            W: array of precomputed twiddle factors
//
// Returns:   values in array x are replaced with result
//
// Calls:     Nothing
//
// Notes:     Bit-reversed address reordering of the sequence
//            is performed in this function.
//            In place computation of FFT. i.e. same memory address
//            that stored the discrete sequence, now contains the FFT
//            output bins.
///////////////////////////////////////////////////////////////////////
{
    cfloat u, temp, tm;
    cfloat *Wptr;

    int i, j, k, len, Windex;

    // perform FFT butterfly
    Windex = 1;
    for(len = n/2 ; len > 0 ; len /= 2) {  // This loop determines the stage of the FFT
        Wptr = W;  // Re-Initializes pointer to start of W
        for (j = 0 ; j < len ; j++) {   // This loop determines which Twiddle factor to use
            u = *Wptr;
            for (i = j ; i < n ; i = i + 2*len) {  // This loops contains the butterfly operation
                                                   // for a particular Twiddle factor in a particular stage
                temp.re = x[i].re + x[i+len].re;
                temp.im = x[i].im + x[i+len].im;
                tm.re = x[i].re - x[i+len].re;
                tm.im = x[i].im - x[i+len].im;

                // The next statements perform complex Multiplication (MAC)
                // (a + ib)(x + iy) = (ax - by) + i(ay + bx)
                x[i+len].re = tm.re*u.re - tm.im*u.im;  // (ax - by)
                x[i+len].im = tm.re*u.im + tm.im*u.re;  // i(ay + bx)
                x[i] = temp;
            }
            Wptr = Wptr + Windex;
        }
        Windex = 2*Windex;
    }

    // Rearrange data by bit reversed addressing
    // this step must occur after the FFT butterfly
    j = 0;
    for (i = 1; i < (n-1); i++) {
        k = n/2;
        while(k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
        if (i < j) {
            temp = x[j];
            x[j] = x[i];   // Here is where the swap happens
            x[i] = temp;
        }
    }  // End of bit-reversed addressing

}  // end of fft_c function

void init_W(int N, cfloat *W)
///////////////////////////////////////////////////////////////////////
// Purpose:   Calculate the twiddle factors needed by the FFT.
//
// Input:     N: length of FFT, W: array to store twiddle factors
//
// Returns:   values are stored in array W
//
// Calls:     Nothing
//
// Notes:     Floats used rather than double to save memory.
///////////////////////////////////////////////////////////////////////
{
    int i;

    float a = 2.0*PI/N;

    for(i = 0 ; i < N ; i++) {
        W[i].re = (float) cos(-i*a);
        W[i].im = (float) sin(-i*a);
    }
}  // end of init_W function
