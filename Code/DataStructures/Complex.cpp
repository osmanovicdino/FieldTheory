#ifndef COMPLEX_CPP
#define COMPLEX_CPP

bool operator>(const complex<double> &a1, const complex<double> &a2) {
    return a1.real() > a2.real();
}

bool operator<(const complex<double> &a1, const complex<double> &a2)
{
    return a1.real() < a2.real();
}

/* 
fftw_complex operator+(const fftw_complex &a1, const fftw_complex &a2) {
    fftw_complex a;
    a[0] = a1[0]+a2[0];
    a[1] = a1[1]+a2[1];

}



fftw_complex operator*(const fftw_complex &a1, double c) {
    fftw_complex b;
    b[0] = a1[0] * c;
    b[1] = a1[1] * c;
    return b;
}

fftw_complex operator*(double c, const fftw_complex &a1)
{
    return c*a1;
}

fftw_complex operator*(const fftw_complex &a, const fftw_complex &b) {
    fftw_complex c;
    c[0] = a[0] * b[0] - a[1] * b[1];
    c[1] = a[1] * b[0] + a[0] * b[1];
    return c;
}
 */
#endif /* COMPLEX_CPP */
