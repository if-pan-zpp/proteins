#include "ran2.hpp"

Scalar util::ran2(Scalar iseed) {
    const int IM1 = 2147483563, IM2 = 2147483399,IMM1 = IM1-1, IA1 = 40014, IA2 = 40692,
        IQ1 = 53668, IQ2 = 52774, IR1 = 12211, IR2 = 3791, NTAB = 32, NDIV = 1+IMM1/NTAB;
    const Scalar EPS = 1.2e-7, RNMX = 1 - EPS, AM = 1 / IM1;

    static int iy, idum2 = 123456789, idum;
    static int iv[NTAB];
    int k, j;

    if(iseed != 0)
        idum = -iseed;
    if(idum <= 0) {
        idum = max(-idum, 1);
        idum2 = idum;
        for(j = NTAB + 8; j>=1; j--) {
            k = idum / IQ1;
            idum = IA1 * (idum - k * IQ1) - k * IR1;
            if(idum < 0)
                idum += IM1;
            if(j - 1 <= NTAB)
                iv[j - 1] = idum;
        }
        iy = iv[0];
    }
    k = idum / IQ1;
    idum = IA1 * (idum - k * IQ1) - k* IR1;
    if(idum < 0)
        idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if(idum2 < 0)
        idum2 += IM2;
    j = 1 + iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;
    if(iy < 1)
        iy += IMM1;
    return min(AM * iy, RNMX);
}