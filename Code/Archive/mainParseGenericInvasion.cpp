#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <random>
//#include <mutex>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <chrono>
#include <type_traits>
//#include <stdafx>
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "DataStructures/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
#include "CahnHilliard/cahnhilliard.h"

#include "fftw3.h"

using namespace std;



int main(int argc, char **argv)
{
    srand(time(NULL));

    double c0 = 0.2;
    double c1 = 0.8;
    double eps = 0.4;




    //line 1 = all parameters;

    //line 2-5 = no. of terms -> rate, pow1, pow2, pow3, pow4 etc
    string importstring;
    string field1s;
    string field2s;
    string field3s;

    double ani;


    if(argc == 6) {
        stringstream ss;
        ss << argv[1];
        importstring = ss.str();

        stringstream ssg1,ssg2,ssg3;
        ssg1 << argv[2];
        ssg2 << argv[3];
        ssg3 << argv[4];
        ani = atof(argv[5]);

        field1s = ssg1.str();
        field2s = ssg2.str();
        field3s = ssg3.str();
    }
    else{
        error("no");
    }

    double T;
    bool err1;
    matrix<double> mat1 = importcsv(importstring,T,err1);

    int n = mat1(0,0);

    vector1<int> phaseseps(n);

    for(int i = 0  ; i < n ; i++)
    phaseseps[i] = mat1(1,i);

    vector1<double> epsi((n)*(n-1)/2);




    for(int i = 0  ; i < (n)*(n-1)/2 ; i++) {
        epsi[i] = mat1(2,i);
    }

    

    int k  = 0;
    matrix<double> epsa(n,n);
    for(int i = 0  ; i < n ; i++) {
        for(int j =  i+1 ; j < n ; j++) {
            epsa(i,j) =  epsi[k];
            epsa(j,i) = epsi[k];
            k++;
        }
    }



    vector1<double> init(n);

    for (int i = 0; i < n; i++)
        init[i] = mat1(3, i);

    // epsilon12 = -mat1(0,4);
    // epsilon13 = -mat1(0,5);
    // epsilon14 = -mat1(0,6);
    // // epsilon12 = mat1(0, 4);
    // // epsilon13 = mat1(0, 5);
    // // epsilon14 = mat1(0, 6);
    // epsilon23 = mat1(0,7);
    // epsilon24 = mat1(0,8);
    // epsilon34 = mat1(0,9);

    // x1 = mat1(0,10);
    // x2 = mat1(0,11);
    // x3 = mat1(0,12);
    // x4 = mat1(0,13);



    // if (argc == 5)
    // {
    //     c0 = atof(argv[1]);
    //     c1 = atof(argv[2]);
    //     eps = atof(argv[3]);
    //     dens = atof(argv[4]);
    // }
    // else {
    //     c0 = 0.2;
    //     c1 = 0.8;
    //     eps = 0.2;
    //     dens = 0.5;
    // }


    typedef complex<double> myc;
    typedef InvasionLinearReversibleB<myc> ILB;
    typedef InvasionLinearReversibleA<myc> ILA;
    typedef Field_Wrapper<myc,myc> FWCC; 
    typedef CoupledPhaseSeparatingSystem<myc> CPSS;
    typedef Rule_Wrapper<myc, myc, myc, myc> RWC;

    CH_builder p;
    int nof = n;
    p.number_of_fields = nof;
    p.N1 = 1024;
    p.N2 = 1024;
    double *str = new double[4];
    str[0] = 0.*1.;
    str[1] = 0.*1.;
    str[2] = 0.*1.;
    str[3] = 0.*1.;
    //GenNoise<myc> b(p);
    modelABnoise A(1.0,1.0);

    CHWithNoise<modelABnoise> a(p,A,str);


    // ILB c0(ratef1, 0.0 * ratef1 /*no growth from below as this is the lowest field*/, 0.0 * rateb1 /*no decay down*/, rateb1, 0, 0, 1, 5);
    // ILB c1(ratef1, ratef1, rateb1, rateb1, 1, 0, 2, 5);
    // ILB c2(ratef1, ratef1, rateb1, rateb1, 2, 1, 3, 5);
    // ILB c3(ratef1, ratef1, rateb1, rateb1, 3, 2, 4, 5);
    // ILB c4(0.0 * ratef1 /*no loss to above*/, ratef1, rateb1, 0.0 * rateb1 /*no decay above*/, 4, 3, 4, 5);


    FWCC my_chemsitry(p);
    NoWeight<myc, myc> nw;
    for (int j = 4 ; j < mat1.getnrows() ; j++)
    {
        int no_chem = mat1(j, 0);
        // cout << no_chem << endl;

        if(no_chem == 0 ) {
            my_chemsitry.add_method(nw,j-4);
        }
        else{
        MultipleReactions<myc> c6(no_chem);

        double tot = 0.0;

        for(int i = 0 ; i < no_chem ; i++) {
            // cout << i << endl;
            vector1<int> jpow(nof);
            for(int k =  i*(nof+1)+2 ; k < i*(nof+1)+2+nof ; k++) {
 
                jpow[k- (i*(nof+1)+2)] = (int)mat1(j,k);
            }

            GenericChemistry<myc> c6_0(mat1(j,i*(nof+1)+1), jpow);
            c6.add_chemical_reaction(c6_0,i);
            double tot1 = 1.0;
            for(int k = 0 ; k < nof ; k++) {
                tot1 *= Power(init[k], jpow[k]);
               // cout << tot1 << ",";
            }
            //cout << endl;

            
            tot += mat1(j, i * (nof + 1) + 1) * tot1;

            // cout << endl;
        }
        //cout << tot << endl;
        
        my_chemsitry.add_method(c6, j-4);
        }
    }

    

    cout << "chemistries added" << endl;
   
    // the solvent;
    // double rate2 = 0.192;
    // double rate1 = 0.712;
    // InvasionLinearReversibleA<myc> c6(-rate1,-rate2,3,2,0);
    // InvasionLinearReversibleA<myc> c7(rate1, rate2, 3, 2, 0);
    // InvasionLinearReversibleA<myc> c8(rate1, rate2, 3, 2, 0);

    
    // double kr1 = 0.514;
    // double kr2 = 0.676;
    // double kf1 = 0.4356;
    // double kf2 = 0.367;

    // MultipleReactions<myc> c6(4);
    // // InvasionSquared<myc> c6_1(kf2,1,0);
    // // InvasionSquared<myc> c6_2(2*kf1, 0, 3);
    // InvasionLinear<myc> c6_1(kf1,0,0);
    // InvasionDecay<myc> c6_2(kf2,0);
    // InvasionLinear<myc> c6_3(-kr1, 0, 1);
    // InvasionLinear<myc> c6_4(-kr2, 2, 3);

    // c6.add_chemical_reaction(c6_1, 0);
    // c6.add_chemical_reaction(c6_2, 1);
    // c6.add_chemical_reaction(c6_3, 2);
    // c6.add_chemical_reaction(c6_4, 3);

    // my_chemsitry.add_method(c6, 0);

    // MultipleReactions<myc> c7(2);
    // InvasionLinear<myc> c7_1(-kf1, 0, 0);
    // InvasionLinear<myc> c7_2(kr1, 0, 1);

    // c7.add_chemical_reaction(c7_1, 0);
    // c7.add_chemical_reaction(c7_2, 1);


    // my_chemsitry.add_method(c7, 1);
    // cout << "done" << endl;
    // MultipleReactions<myc> c8(2);
    // InvasionDecay<myc> c8_1(-kf2, 0);
    // InvasionLinear<myc> c8_2(kr2, 2, 3);

    // c8.add_chemical_reaction(c8_1, 0);
    // c8.add_chemical_reaction(c8_2, 1);

    // my_chemsitry.add_method(c8, 2);
    // cout << "done" << endl;
    // MultipleReactions<myc> c9(2);
    // InvasionDecay<myc> c9_1(-kf2, 0);
    // InvasionLinear<myc> c9_4(kr2, 2, 3);

    // c9.add_chemical_reaction(c9_1, 0);
    // c9.add_chemical_reaction(c9_4, 1);

    // my_chemsitry.add_method(c9, 3);

    // cout << "added chemistries" << endl;

    // myc ** calculated_reactions = new myc *[p.number_of_fields];

    // for(int k =0 ; k < p.number_of_fields ; k++) {
    // calculated_reactions[k] = (myc *)fftw_malloc(p.N1 * p.N2 * sizeof(myc));
    // for (int j = 0; j < p.N1 * p.N2; j++)
    // {
    //     calculated_reactions[k][j] = 1.0;
    // }
    // }
    // my_chemsitry.Calculate_Results(calculated_reactions);

    // cout << "up to here" << endl;

    FWCC my_weights(p);

    double nu = 1.0;


    {
        int k = 0;
        vector1<int> which1(nof-1);
        vector1<double> myeps(nof-1);
        for(int i = 0  ;  i < nof ; i++) {
            if(i==0) {

            }
            else{
                which1[k] = i;
                myeps(k) = epsa(0,i);
                k++;
            }
        }
        CahnHilliardWithCouplingWeightGenericNSQR<myc> w0(c0,c1, nu,myeps,which1);
        my_weights.add_method(w0, 0);
    }


    for(int lk = 1 ; lk < nof ; lk++) {
        int k = 0;
        vector1<int> which1(nof - 1);
        vector1<double> myeps(nof - 1);
        for (int i = 0; i < nof; i++)
        {
            if (i == lk)
            {
            }
            else
            {
                which1[k] = i;
                myeps(k) = epsa(lk, i);
                k++;
            }
        }

        DiffDiffusiveWeightGenericNSQR<myc> w1(myeps,which1);

        my_weights.add_method(w1, lk);
    }


    //my_weights.add_method(w0, 0);

    // my_weights.add_method(w1, 1);

    // my_weights.add_method(w2, 2);

    // my_weights.add_method(w3, 3);

    cout << "weights added" << endl;

    a.set_chems(my_chemsitry);
    a.set_weights(my_weights);




    double dt = 0.05;
    //double dx = 0.05;

    //double temp1 = (1. / (dx * p.N1)) ;

    double surf = eps;
    double L = 100.0;
    double temp1 = SQR(2.*pii/L);

    // //cout << (1. / (dx * p.N1)) << endl;
    // cout << temp1 << endl;

    double diffusion_constant =  1.;


    // //double surface_width = 2.0;

    RWC my_rules(p);
    // //MultiPhaseBundleComplex Bund(w0, p, dt, Di, temp1, surface_width, surf);

    DiffusionWithSurfaceTension<myc> C(p, dt, diffusion_constant, temp1, surf);
    my_rules.add_method(C, 0); // as this is a collective one, only need to define the first rule;


    for(int lk = 1 ; lk < nof ; lk++) {
        double diff = diffusion_constant;
        if(lk>=3) {
             diff *=100.;
        }
    DiffusionWithInteraction<myc> C1(p, dt, diff, temp1);
    my_rules.add_method(C1, lk);
    }
    // DiffusionWithInteraction<myc> C2(p, dt, diffusion_constant, temp1);
    // DiffusionWithInteraction<myc> C3(p, dt, diffusion_constant, temp1);
    // // MultiPhase r1(p, 6, dt, Di, temp1, surface_width, ChiMatrix);
    // // MultiPhase r2(p);
    // // MultiPhase r3(p);
    // // MultiPhase r4(p);
    // // MultiPhase r5(p);
    // // MultiPhase r6(p); //all empty

    // // cout << "created diffusion" << endl;

    // //for(int k = 0  ; k < nof ; k++)

    // my_rules.add_method(C1, 1); //as this is a collective one, only need to define the first rule;
    // my_rules.add_method(C2, 2); //as this is a collective one, only need to define the first rule;
    // my_rules.add_method(C3, 3); //as this is a collective one, only need to define the first rule;

    cout << "all rules added" << endl;
    // cout << "all methods added" << endl;


    a.set_rules(my_rules);

    // cout << "Setting rules" << endl;

    cout << "old rules set" << endl;



    // for(int k = 0  ; k < nof ; k++)
    //     a.rules.chems[k]->print();

    // for (int k = 0; k < nof; k++)
    //     a.newrules.chems[k]->print();

    // cout << "added methods" << endl;

    vector<matrix<myc> > v;

    for(int j = 0 ; j < nof ; j++) {
    matrix<myc> field1(p.N1, p.N2);
    v.push_back(field1);
    }

    matrix<double> f1 = importcsv(field1s, T, err1);
    matrix<double> f2 = importcsv(field2s, T, err1);
    matrix<double> f3 = importcsv(field3s, T, err1);
    // double x1 = (kr2*sqrt((Power(kf2,2)*kr1*Power(x2,5))/(kf1*Power(kr2,2)*Power(x4,4)))*Power(x4,2))/(kf2*Power(x2,2));
    // double x3 = sqrt((Power(kf2, 2) * kr1 * Power(x2, 5)) / (kf1 * Power(kr2, 2) * Power(x4, 4)));
    double maxval;
    (f1+f2).maxima(maxval);
    double meanofinv = 0.0;
    // for(int lk  = 0  ; lk < nof ; lk++) {
    //     double x1 = init[lk];
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                //double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[0](i, j) = f1(i,j);
                v[1](i, j) = f2(i, j);
                v[2](i, j) = f3(i, j);
                v[3](i, j) = 1.E-6 + maxval - (f1(i, j))-f2(i,j);
                meanofinv += v[3](i,j).real();
                v[4](i, j) =  0.0;
            }
        }
    //}

        meanofinv /= SQR((double)(p.N1));

        // cout << meanofinv << endl;
        myc fac = (ani / meanofinv);
        v[3] *= fac;

        // double val;
        // (field1+(5.*field2)+field3).maxima(val);
        // if(val>1.) error("initial packing fracs not correct");

        for (int lk = 0; lk < nof; lk++)
        {
            a.set_field(v[lk], lk);
    }



    // my_chemsitry.Calculate_Results(a.fields);

    // for(int i = 0  ; i < nof ; i++)
    // cout << a.fields[i][0] <<",",
    // cout << endl;

    // for (int i = 0; i < nof; i++)
    //     cout << my_chemsitry.calculated_reactions[i][0] << ",",
    //         cout << endl;
    // pausel();
    // for(int k = 0  ; k < nof ; k++) {
    // cout << a.fields[k][0] << ",";
    // }
    // cout << endl;

    // for (int k = 0; k < nof; k++)
    // {
    //     cout << a.fields[k][1000] << ",";
    // }
    // cout << endl;

    // for (int i = 0; i < nof; i++)
    //     cout << a.fields[i][0] << endl;
    // a.chems.Calculate_Results(a.fields);
    // for(int i = 0  ; i < nof ; i++)
    //     cout << a.fields[i][0] << endl;
    // cout << "chem done" << endl;
    // pausel();

    // a.set_field(field2, 1);

    // a.set_field(field3, 2);

    // a.set_field(field4, 3);

    cout << "all fields set" << endl;
 

    cout << init << endl;
    // auto start = std::chrono::high_resolution_clock::now();
    int runtime = 10001;
    int every = 100;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    // auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < runtime; i++)
    {

        if (i % every == 0 && i > 0)
        {
            // stringstream strep1;
            // stringstream strep2;
            // stringstream strep3;
            // stringstream strep4;

            // strep1 << dens;
            // strep4 << c0;
            // strep2 << c1;
            // strep3 <<  surf;
            string s1 = importstring;
            // string s1 = "denp=" + strep1.str() + "c0=" + strep4.str() + "_c1=" + strep2.str() + "_surf=" + strep3.str();
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << i / every;
            string s2 = "_i=" + ss.str();
            string su = s1.substr(0, s1.size() - 4);
            cout << su + s2 << endl;
            a.print_all_results(su + s2);
        }
        cout << i << endl;
        a.Update();
    }



    // auto stop = std::chrono::high_resolution_clock::now();

    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // cout << duration.count() << endl;

    // cout << "he's done ya again" << endl;

    /*
    int N1 = 512;
    int N2 = 512;
    double **an_in_array = new double *[4];
    double **an_out_array = new double *[4];
    an_in_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));
    an_out_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));

    // fftw_plan p;

    auto start = std::chrono::high_resolution_clock::now();

    // p = fftw_plan_r2r_2d(N1, N2, an_in_array[0], an_out_array[0], FFTW_REDFT10, FFTW_REDFT10, 1);

    for(int i = 0 ; i < 1000 ; i++) {
        fftw_plan p;

        p = fftw_plan_r2r_2d(N1, N2, an_in_array[0], an_out_array[0], FFTW_REDFT10, FFTW_REDFT10, 1);

        fftw_execute(p);

        fftw_destroy_plan(p);
    }

    auto stop =  std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

    cout << duration.count() << endl;
    */

    // fftw_destroy_plan(p);

    /*
CH_builder p;
p.number_of_fields = 4;
p.N1=1024;
p.N2=1024;



Field_Wrapper a(p);


InvasionAntiInvasionLinear b(1.0,0.5);
NoWeight c;

a.add_method(b,0);
a.add_method(c, 1);
a.add_method(c, 2);
a.add_method(c, 3);

double **fields = new double *[4];
fields[0] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));
fields[1] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));
fields[2] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));
fields[3] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));

for(int i = 0 ; i < 4 ; i++) {
    for(int j = 0 ; j < SQR(1024) ; j++) {
        fields[i][j] = 0.2;
    }
}

cout << "fine to here" << endl;

a.Calculate_Results(fields);

cout << a.calculated_reactions[0][0] << endl;
cout << "done" << endl;
pausel();
*/
    //Reaction_Wrapper()

    // CH a(p);

    // a.chems[0]->what_am_I();
    // a.update_rules[0]->what_am_I();

    // InvasionAntiInvasion *b = new InvasionAntiInvasion;

    // cout << a.chems.size() << endl;

    // a.chems[0] =  b;

    // a.chems[0]->what_am_I();

    /*
    int N1 = 512;
    int N2 = 512;
    double **an_in_array = new double*[4];
    double **an_out_array = new double*[4];
    an_in_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));
    an_out_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));

    fftw_plan p;

    p = fftw_plan_r2r_2d(N1, N2, an_in_array[0], an_out_array[0], FFTW_REDFT10, FFTW_REDFT10, 1);

    fftw_execute(p);
    */

    /*
int N1 = 512;
int N2 = 512;



double *an_array;
an_array = (double *)fftw_malloc(N1 * N2 * sizeof(double));

double *out;
out = (double *)fftw_malloc(N1 * N2 * sizeof(double));
// in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
// out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

for (int i = 0; i < N1; i++)
{
    for (int j = 0; j < N2; j++)
    {
        an_array[i * N1 + j] = cos(pi*i*j/10.);
    }
}

fftw_plan p;

p = fftw_plan_r2r_2d(N1, N2, an_array, out, FFTW_REDFT10, FFTW_REDFT10, 1);

fftw_execute(p);
cout << an_array[0] << endl;
cout << out[0] / (4 * 512) << endl;




outfunc2D(out,N1,N2,"cosine");

fftw_plan p2;
p2 = fftw_plan_r2r_2d(N1, N2, out, an_array, FFTW_REDFT01, FFTW_REDFT01, 1);

fftw_execute(p2);

cout << an_array[0]/SQR(2*512) << endl;
cout << an_array[1] << endl;
cout << out[0] / (4*512) << endl;
*/

    // cout << a <<endl;
}