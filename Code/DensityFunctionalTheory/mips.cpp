#ifndef MIPS_CPP
#define MIPS_CPP



void outfunc(vecmat &a, string s) {
    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    int ni  = a.getsize();
    int nj  = a[0].getnrows();
    int nk  = a[0].getncols();

    for(int i = 0 ; i < ni ; i++) {
        for(int j = 0  ; j < nj ; j++) {
            for(int k = 0  ; k < nk ; k++) {
                if(j == nj-1 && k == nk-1) {
                    myfileg << a[i](j,k);
                }
                else{
                    myfileg << a[i](j, k) << ",";
                }
            }
        }
        myfileg << "\n";
    }

    myfileg.close();
}

mips_mf::mips_mf(int Nrr, int Nqq, int NN, double ll, double v00) : N(NN), Nr(Nrr), Nq(Nqq), dat(vector1<matrix<double>>(Nqq)), potm(vector1<matrix<COMPLEX>>(Nq * Nq)), l(ll), xpos(vector1<double>(Nrr)), ypos(vector1<double>(Nrr)), qpos(vector1<double>(Nqq))
{

dx = l/(double)(Nrr);

for(int i = 0 ; i < Nrr ; i++) {
    xpos[i] = i*dx;
    ypos[i] = i*dx;
}
dq = (2.*pi)/(double)(Nqq);

for(int i =0 ; i < Nqq; i++) {
    qpos[i] = i*dq;
}


double dk  =  (2*pi/l);


for(int i = 0  ; i < Nqq ; i++) {
    matrix<double> res(Nrr,Nrr,(double)N/(l*l*2*pi));
    dat[i]=res;
}

// for(int i = 0  ; i < Nr*Nr ; i++) {
//     matrix<double> res(Nrr, Nrr);
//     int i1 = i/Nr;
//     int i2 = i%Nr;

//     double kx =  dk*i1;
//     double ky  = dk*i2;

//     for(int j = 0  ; j < Nr ; j++) {
//         for(int k = 0 ; k < Nr ; k++) {
//             res(j,k) = (1./(double)Nr)* cos(kx*xpos[j] + ky*ypos[k]);
//         }
//     }

    
//     fm[i] = res;
// }

epsilon = 1.0;
sigma = 1.0;
v0 = v00;


double x0 =  xpos[0];
double y0 = ypos[0];

vecmat temp(Nq);

for(int i = 0 ; i < Nq*Nq ; i++) {
    //matrix<double> res(Nrr, Nrr);
    int q1 =  i/Nq;
    int q2 =  i%Nq;

    matrix<COMPLEX> potmat(Nr,Nr);

    double tot = 0.0;
    for(int j = 0  ; j < Nr ; j++ ) {
        for(int k = 0 ; k < Nr ; k++) {
            double hm = potential(xpos[j], ypos[k], qpos[q1], x0, y0, qpos[q2]);
            potmat(j,k).real = hm;
            tot += hm;
        }
    }


    
    if (Nr != 0 && (Nr & (Nr - 1)) == 0) {
        FFT2D(potmat, Nr, Nr, 1);
    }

    potm[i] = potmat;
}


}

// matrix<double> mips_mf::fourier(matrix<double> &m1) {
//     matrix<double> res(Nr,Nr);
//     #pragma omp prallel for
//     for(int kx = 0 ; kx < Nr ; kx++) {
//         for(int ky  = 0 ; ky < Nr ; ky++) {
//             res(kx,ky)=sum_all_elements(fm[kx*Nr+ky]&m1);
//         }
//     }
//     return res;
// }

double mips_mf::min(vector1<matrix<double > > &a) {
    vector1<double> minp(a.getsize());
    for(int i = 0  ; i < a.getsize() ; i++) {
        double b;
        a[i].minima(b);
        minp[i] =  b;
    }
    return minval(minp);
}

double mips_mf::max(vector1<matrix<double>> &a)
{
    vector1<double> minp(a.getsize());
    for (int i = 0; i < a.getsize(); i++)
    {
        double b;
        a[i].maxima(b);
        minp[i] = b;
    }
    return maxval(minp);
}

// inline double mips_mf::potential(double x1, double y1, double z1, double x2, double y2, double z2)
// {

//     double dx = x1 - x2;
//     double dy = y1 - y2;
//     if (abs(dx) > l / 2.)
//         dx = dx - SIGN(l, dx);

//     if (abs(dy) > l / 2.)
//         dy = dy - SIGN(l, dy);

//     double r2 = SQR(dx) + SQR(dy);

//     double fac1 = (sigma / sqrt(r2));
//     double fac2 = SQR(fac1);
//     //else expf = 0.0;

//     double st;
//     if (r2 > 0.69 * SQR(sigma))
//     { //cut off
//         double f1 = ((4. * epsilon * v0)) * (4 * fac2 * fac2 - 2 * fac2); //long range potential
//         st = f1 * (dx * cos(z1) - dx * cos(z2) + dy * sin(z1) - dy * sin(z2));
//     }
//     else
//     {
//         st = 0.0;
//     }

//     return st;
// }

inline double mips_mf::potential(double x1, double y1, double z1, double x2, double y2, double z2)
{

    double dx = x1 - x2;
    double dy = y1 - y2;
    if (abs(dx) > l / 2.)
        dx = dx - SIGN(l, dx);

    if (abs(dy) > l / 2.)
        dy = dy - SIGN(l, dy);

    double r2 = SQR(dx)+SQR(dy);                    //cut off
    double f1 = ((4. * epsilon * v0)) * exp(-r2); //long range potential
    double st = f1 * (dx * cos(z1) - dx * cos(z2) + dy * sin(z1) - dy * sin(z2));


    return st;
}


// inline double mips_mf::potential(double x1, double y1, double z1, double x2, double y2, double z2)
// {

//     double dx = x1 - x2;
//     double dy = y1 - y2;
//     if (abs(dx) > l / 2.)
//         dx = dx - SIGN(l, dx);

//     if (abs(dy) > l / 2.)
//         dy = dy - SIGN(l, dy);

//     double r2 = SQR(dx) + SQR(dy);


//     return epsilon*exp(-r2/sigma);
// }

// void mips_mf::gen_fm(double kx, double ky, matrix<double> &m1, matrix<double> &m2) {
//     for(int i = 0  ; i < Nr ; i++)
//         for(int j = 0  ; j < Nr ; j++) {
//             m1(i, j) = cos(kx * xpos[i] + ky * ypos[j]);
//             m2(i, j) = sin(kx * xpos[i] + ky * ypos[j]);
//         }
// }

// void mips_mf::fourier() {

//     for(int kx = 0 ; kx < Nr ; kx++) {
//         for(int ky  = 0 ; ky < Nr ; ky++) {
//             matrix<double> m1(Nr,Nr);
//             matrix<double> m2(Nr, Nr);
//             gen_fm(kx,ky)
//         }
//     }
// }
matrix<double> mips_mf::qintegral()
{
    matrix<double> fd(Nr, Nr);
    for (int i = 0; i < Nq; i++)
        fd += dat[i];
    return dq*fd;
}

double mips_mf::integrate()
{
    double tot  = 0.0;
    matrix<double> fd(Nr, Nr);
    for (int i = 0; i < Nq; i++)
        for(int j = 0  ; j < Nr ; j++)
            for(int k = 0; k < Nr ; k++)
                tot += dat[i](j,k);
    return dq * dx * dx * tot;
}

double mips_mf::integrate(vecmat &a)
{
    double tot = 0.0;
    matrix<double> fd(Nr, Nr);
    for (int i = 0; i < Nq; i++)
        for (int j = 0; j < Nr; j++)
            for (int k = 0; k < Nr; k++)
                tot += a[i](j, k);
    return dq * dx * dx * tot;
}

matrix<double> mips_mf::calc_HS() {

matrix<double> fd(Nr,Nr);

for(int i = 0  ; i < Nq ; i++)
    fd += dq*dat[i];

fd *= (pi/4.)*SQR(sigma);


matrix<double> res(Nr,Nr);

for(int i = 0  ; i < Nr ; i++)
for(int j = 0 ; j < Nr ; j++) {
    res(i,j) = -log(1-fd(i,j))+(3*fd(i,j)-2*SQR(fd(i,j)))/SQR(1-fd(i,j));
}

return res;

}


vector1< matrix<double> > mips_mf::calc_int()
{

    vector1<matrix<double> > inter(Nq);

    // \int dr rho(,x,y,z) phi(x,y,z,x1,y1,z1)
    int maxr = 10;
    #pragma omp parallel for
    for(int i = 0 ; i < Nq ; i++) {
        inter[i] = matrix<double>(Nr,Nr);
        for(int s1 = 0 ; s1 < Nr ; s1++) {
            for(int s2 = 0 ; s2 < Nr ; s2++) {
                //cout << i << " " << s1 << " " << s2 << endl;
                double sum = 0.0;
                double x1 = xpos[s1];
                double y1 = ypos[s2];
                double q1 = qpos[i];

                for (int i2 = 0; i2 < Nq; i2++)
                {
                    for (int s12 = s1-maxr; s12 < s1+maxr; s12++)
                    {
                        for (int s22 = s2 - maxr; s22 < s2+maxr; s22++)
                        {
                            int s12t = s12;
                            int s22t = s22;
                            if(s12t < 0) s12t += Nr;
                            else if(s12t > Nr-1) s12t -= Nr;
                            else {}

                            if (s22t < 0)
                                s22t += Nr;
                            else if (s22t > Nr - 1)
                                s22t -= Nr;
                            else{}

                            sum += dx*dx*dq*dat[i2](s12t,s22t) * potential(x1,y1,q1,xpos[s12t],ypos[s22t],qpos[i2]);
                        }
                    }
                }
                inter[i](s1,s2) = sum;
            }
        }
    }

    return inter;
    
}

vector1<matrix<double>> mips_mf::calc_int2()
{

    vector1<matrix<double>> inter(Nq);
    #pragma omp parallel for
    for(int i = 0  ; i < Nq ; i++) {
        inter[i] = matrix<double>(Nr,Nr);
        matrix<COMPLEX> fin(Nr,Nr);

        for(int j = 0 ; j < Nq ; j++) {
            matrix<COMPLEX> den2(dat[j]);
            FFT2D(den2, Nr,Nr, 1);
            matrix<COMPLEX> res = den2&potm[i*Nq+j];
            FFT2D(res,Nr,Nr,-1);

            fin += res;

            
        }
        for(int k =  0 ; k < Nr ; k++)
            for(int h = 0 ; h < Nr ; h++)
                //inter[i](k, h) = (dq * dx *  dx* Nr * Nr) *fin(k, h).getreal();
                inter[i](k, h) = (dq * Nr * Nr) * fin(k, h).getreal();
    }
    return inter;

    // \int dr rho(,x,y,z) phi(x,y,z,x1,y1,z1)
  
}

void mips_mf::updatedat(double dt) {
    
    matrix<double> hs = this->calc_HS();
    vector1<matrix<double> > phi = this->calc_int2();

    double norm  = 0.0;
    double minv,maxv;
    hs.minima(minv);
    hs.maxima(maxv);

    cout << "min: " << this->min(dat) <<endl;
    cout << "max: " << this->max(dat) << endl;

    cout << "minf: " << this->min(phi) << endl;
    cout << "maxf: " << this->max(phi) << endl;

    //  cout << hs << endl;
    cout << "minHS: " << minv << endl;
    cout << "maxHS: " << maxv << endl;

    outfunc(this->qintegral(), "int");
    outfunc(dat,"den");
    outfunc(phi, "phi");


    for (int i = 0; i < Nq; i++)
    {
        
        for (int s1 = 0; s1 < Nr; s1++)
        {
            for (int s2 = 0; s2 < Nr; s2++)
            {
                double upd = log((dat[i](s1,s2)) )+ ((1./(2.*pi))*hs(s1,s2)+phi[i](s1,s2)) ;
               
                dat[i](s1,s2) = dat[i](s1,s2) - dt*upd;
                norm += dat[i](s1, s2) * dx * dx * dq;
            }
        }
    }

    dat *= ((double)N / norm);
}

void mips_mf::updatedatstore()
{

    matrix<double> hs = this->calc_HS();
    vector1<matrix<double>> phi = this->calc_int2();

    double norm = 0.0;
    double minv, maxv;
    hs.minima(minv);
    hs.maxima(maxv);

    cout << "min: " << this->min(dat) << endl;
    cout << "max: " << this->max(dat) << endl;

    cout << "minf: " << this->min(phi) << endl;
    cout << "maxf: " << this->max(phi) << endl;

    //  cout << hs << endl;
    cout << "minHS: " << minv << endl;
    cout << "maxHS: " << maxv << endl;


    vecmat temp(Nq);


    for (int i = 0; i < Nq; i++)
    {
        temp[i]=matrix<double>(Nr,Nr);
        for (int s1 = 0; s1 < Nr; s1++)
        {
            for (int s2 = 0; s2 < Nr; s2++)
            {
                double upd = log(dat[i](s1, s2))+((1. / (2. * pi))  * hs(s1, s2) + phi[i](s1, s2));
                temp[i](s1,s2) = upd;
                // dat[i](s1, s2) = dat[i](s1, s2) - dt * upd;
                // norm += dat[i](s1, s2) * dx * dx * dq;
            }
        }
    }

    outfunc(this->qintegral(), "int");
    outfunc(dat, "den");
    outfunc(phi, "phi");
    outfunc(temp, "fun");

    double fac = this->integrate(temp);


    double lambda = fac/(2*pi*l*l);

    //cout << fac-lambda*(this->integrate()) << endl;

    double min = this->min(temp) - lambda;
    double max = this->max(temp) - lambda;

    double minp = this->min(dat);
    double maxp = this->max(dat);

    // cout << maxp << endl;
    // cout << minp << endl;
    // cout << min << endl;
    // cout << max << endl;
    // cout << lambda << endl;

    double dt1 = minp/max;
    double dt2 = ((-4./pi)+(1./(2.*pi))*maxp)/(min);

    double dt = 0.5*(MIN(dt1,dt2));
    //cout << this->integrate() << endl;

    for (int i = 0; i < Nq; i++)
    {
        for (int s1 = 0; s1 < Nr; s1++)
        {
            for (int s2 = 0; s2 < Nr; s2++)
            {
                
                dat[i](s1, s2) = dat[i](s1, s2) - dt * (temp[i](s1, s2) - lambda);

            }
        }
    }
    cout << "dt: " << dt << endl;
    


    //dat *= ((double)N / norm);
}

#endif /* MIPS_CPP */
