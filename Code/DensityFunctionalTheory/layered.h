#ifndef LAYERED_H
#define LAYERED_H

#include "../DataStructures/vector1.h"

#include "../DataStructures/matrix2.h"


struct layered {

int Nx,Ny,Nz;

double dx,dy,dz;

double na,nb,nc;
double chi12,chi13,chi23;
double l12,l13,l23;



int dimension;

vector1<double> density1;
vector1<double> density2;

layered();

double chi(double wab,double waa, double wbb) {
return 3*(wab-0.5*(waa+wbb));
}

void setchis(double w12, double w11, double w22) {
chi12 = chi(w12, w11, w22);
chi13 = chi(0., w11, 0.);
chi23 = chi(0., w22, 0.);


}

void neighboring1D(int i, int &im, int &ip)
{
    if(i==0) {
    im = Nx-1;
    ip = 1;
    }
    else if(i==Nx-1) {
    im = Nx-2;
    ip = 0;
    }
    else{
        im = i-1;
        ip = i+1;
    }
}

void neighboring2D(int i, int j, int &im, int &ip, int &jm, int &jp) {
    int im1,ip1;
    int jm1, jp1;

    neighboring1D(i, im1, ip1);
    neighboring1D(j, jm1, jp1);


    im = im1*Nx + j;
    ip = ip1*Nx + j;

    jm = i*Nx + jm1;
    jp = i *Nx + jp1;
}

void neighboring3D(int i, int j, int k, int &im, int &ip, int &jm, int &jp, int &km, int &kp) {
    int im1, ip1;
    int jm1, jp1;
    int km1, kp1;

    neighboring1D(i, im1, ip1);
    neighboring1D(j, jm1, jp1);
    neighboring1D(k, km1, kp1);

    im = im1 * Nx * Ny + j * Ny + k;
    ip = ip1 * Nx * Ny + j * Ny + k;

    jm = i * Nx * Ny + jm1 * Ny + k;
    jp = i * Nx * Ny + jp1 * Ny + k;

    km = i * Nx * Ny + j * Ny + km1;
    kp = i * Nx * Ny + j * Ny + kp1;
}

void update2D();

void update3D();


};

#include "layered.cpp"

#endif /* MIPS_H */
