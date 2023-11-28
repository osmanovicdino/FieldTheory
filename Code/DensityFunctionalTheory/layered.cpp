#ifndef LAYERED_CPP
#define LAYERED_CPP

layered::layered() : density1(vector1<double>(128 * 128 * 128)), density2(vector1<double>(128 * 128 * 128))
{
    Nx = 128;
    Ny = 128;
    Nz = 128;

    dx = 1.;
    dy = 1.;
    dz = 1.;

    na = 1.0;
    nb = 1.0;
    nc = 1.0;

    chi12 = 0;
    chi13 = 0;
    chi23 = 0;

    l12=1.0;
    l13=1.0;
    l23=1.0;

    dimension = 3;
}

void layered::update3D()
{
    vector1<double> d1(Nx * Ny * Nz);
    vector1<double> d2(Nx * Ny * Nz);
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {

                int im, ip, jm, jp, km, kp;
                neighboring3D(i, j, k, im, ip, jm, jp, km, kp);

                double cxm1 = density1[im];
                double cxp1 = density1[ip];

                double cym1 = density1[jm];
                double cyp1 = density1[jp];

                double czm1 = density1[km];
                double czp1 = density1[kp];

                double cxm2 = density2[im];
                double cxp2 = density2[ip];

                double cym2 = density2[jm];
                double cyp2 = density2[jp];

                double czm2 = density2[km];
                double czp2 = density2[kp];

                double c1 = density1[i * Nx * Ny + j * Nx + k];
                double c2 = density2[i * Nx * Ny + j * Nx + k];

                d1[i * Nx * Ny + j * Nx + k] = (cxm1 + cxp1 - 2 * c1) / (dx * dx) + (cym1 + cyp1 - 2 * c1) / (dy * dy) + (czm1 + czp1 - 2 * c1) / (dz * dz);

                d2[i * Nx * Ny + j * Nx + k] = (cxm2 + cxp2 - 2 * c2) / (dx * dx) + (cym2 + cyp2 - 2 * c2) / (dy * dy) + (czm2 + czp2 - 2 * c2) / (dz * dz);
            }
        }
    }

    vector1<double> upd1(Nx * Ny * Nz);
    vector1<double> upd2(Nx * Ny * Nz);

    for (int i = 0; i < Nx * Ny * Nz; i++)
    {

        double c1 = density1[i];
        double c2 = density2[i];

        upd1[i] = (1. / na) - (1. / nc) + (log(c1) / na) - (log(1 - c1 - c2) / nc) + chi13 - 2 * chi13 * c1 - (chi13 + chi23 - chi12) * c2 - l13 * chi13 * d1[i] - 0.5 * (-l12 * chi12 + l13 * chi13 + l23 * chi23) * d2[i];
        upd2[i] = (1. / nb) - (1. / nc) + (log(c2) / nb) - (log(1 - c1 - c2) / nc) + chi23 - 2 * chi23 * c2 - (chi13 + chi23 - chi12) * c1 - l23 * chi23 * d2[i] - 0.5 * (-l12 * chi12 + l13 * chi13 + l23 * chi23) * d1[i];
    }

    double lambda1 = trap(upd1, dx * dy * dz) / (dx * dy * dz * Nx * Ny * Nz);
    double lambda2 = trap(upd2, dx * dy * dz) / (dx * dy * dz * Nx * Ny * Nz);

    // find the maximally negative value
    double mindt = 10.;
    for (int i = 0; i < Nx * Ny * Nz; i++)
    {
        if (upd1[i] - lambda1 > 0)
        {
            // cout << upd1[i] - lambda1 << endl;
            double dt = density1[i] / (upd1[i] - lambda1);
            if (dt < mindt)
            {

                mindt = dt;
            }
        }

        if (upd2[i] - lambda2 > 0)
        {
            // cout << upd2[i] - lambda2 << endl;
            double dt = density2[i] / (upd2[i] - lambda2);
            if (dt < mindt)
            {

                mindt = dt;
            }
        }

        double dtx = (-1 + density1[i] + density2[i]) / (upd1[i] + upd2[i] - lambda1 - lambda2);

        if(dtx > 0) {
            if (dtx < mindt)
            {

                mindt = dtx;
            }
        }
    }

    

    double dt2 = 0.7 * mindt;
    //cout << dt2 << endl;

    for (int i = 0; i < Nx * Ny * Nz; i++)
    {
        density1[i] = density1[i] - dt2 * (upd1[i] - lambda1);
        density2[i] = density2[i] - dt2 * (upd2[i] - lambda2);
    }
}

#endif /* MIPS_H */
