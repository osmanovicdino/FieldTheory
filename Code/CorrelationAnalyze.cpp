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
#include <algorithm>
#include <parallel/algorithm>
#include <string.h>
#include <dirent.h>

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
// #include "CahnHilliard/cahnhilliard.h"

//#include "fftw3.h"

using namespace std;

bool conv(double a) {
    if(a>0.4) return true;
    else return false;
}

int return_csv_in_dir(string directory, string match, vector<string> &files)
{
    DIR *dir;
    struct dirent *diread;
    // vector<string> files;

    cout << "beginning directory read" << endl;

    if ((dir = opendir(directory.c_str())) != nullptr)
    {
        while ((diread = readdir(dir)) != nullptr)
        {
            std::string fname = diread->d_name;
            if (fname.find(match) != std::string::npos && fname.find(".csv") != std::string::npos)
                files.push_back(fname);

            // files.push_back(diread->d_name);
        }
        closedir(dir);
    }
    else
    {
        perror("opendir");
        return EXIT_FAILURE;
    }

    sort(files.begin(), files.end());

    return 0;
}

int periodicF(int i, int N) {
    if(i>N) return i-N; 
    else return i;
}

double pearson(vector1<double> &v1, vector1<double> &v2) {
    double m1 =  meanish(v1);
    double m2 = meanish(v2);

    double s1=0;
    double s2=0;
    double c1=0;


    for(int i =0 ; i < v1.getsize() ; i++) {
        c1+=(v1[i]-m1)*(v2[i]-m2);
        s1+=SQR(v1[i]-m1);
        s2+=SQR(v2[i]-m2);
    }
    // cout << v1.getsize() << endl;
    // cout << c1 << " " << s1 << " " << s2 << endl;
    // pausel();
    return c1/sqrt(s1*s2);
}

int main(int argc, char **argv)
{
    uint64_t microseconds_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // cout << microseconds_since_epoch << endl;
    // cout << time(NULL) << endl;
    int seed = microseconds_since_epoch % time(NULL);
    srand(seed);

    string mydir;
    string substr;
    string strp;
    int spaceoffset;
    int timeoffset;
    int range1;
    int range2;
    // string mydir = "/u/scratch/d/dinoo/GrowthRun1/";

    if (argc == 4)
    {
        mydir = string(argv[1]);
        substr = string(argv[2]);
        strp = string(argv[3]);
    }
    else
    {
        error("incorrect number of arguments for this executable");
        // cout << "on your own" << endl;
    }

    vector<string> bindfiles;
    cout << mydir << endl;
    return_csv_in_dir(mydir, substr, bindfiles);

    int n = bindfiles.size();
    cout << n << endl;


    // double T;
    // bool err1;
    // bool gen;
    // matrix<bool> a = importcsv("/home/dino/Documents/Chemistry/SubDiffusion/Simulations/LivingDropletsNonNegativeChangeLong/sim5/field0res_chem=7g_i=0001_real.csv",T,gen,err1,conv);

    // cout << a << endl;
    vector<matrix<double> > alldata;
    for(int filen =  0 ; filen < n ; filen++) {
    double T;
    bool err1;
    bool gen;
    cout << bindfiles[filen] << endl;
    matrix<double> bindtemp = importcsv(mydir + "/" + bindfiles[filen],T,err1);
    alldata.push_back(bindtemp);
    }

    // int p00 = 0;
    // int p01 = 0;
    // int p10 = 0;
    // int p11 = 0;



    // int p0x = 0;
    // int p1x = 0;

    // int p0y = 0;
    // int p1y = 0;

    int NR = alldata[0].getnrows();

    int NC = alldata[0].getncols();


    int samp = 1000;
    ofstream myfilex;
    myfilex.open(strp.c_str());
    // matrix<double> mia(spaceoffset+1,timeoffset+1);
    for(int x1 = 0  ; x1 < 200 ; x1++) {
        for(int y1 = 0 ; y1 < 200 ; y1++) {
            for(int t1 = 0 ; t1 < 1; t1+=2) {
                vector1<double> v1(samp);
                vector1<double> v2(samp);
                int iter = 0;
                for(int k = 0 ; k < samp ; k++) {
                    int t2 = (rand() % (n - t1));
                    int x2 = (rand() % (NR - x1));
                    int y2 = (rand() % (NC - y1));
                    //cout << t2 << " " << x2 << " " << y2 << endl;
                    v1[iter] = alldata[t2](x2,y2);
                    v2[iter] = alldata[t2+t1](x2+x1,y2+y1); 
                    iter++;
                }
                myfilex << x1 << "," << y1 << "," << t1 << "," << pearson(v1, v2) << endl;
            }
        }
    }
    myfilex.close();
    

    // for(int timeoffset1 = 0 ; timeoffset1 <= timeoffset ; timeoffset1+=5) {
    // for(int spaceoffset1 = 0 ; spaceoffset1 <= spaceoffset; spaceoffset1+=5) {
    //     matrix<int> pxy(2, 2);
    //     vector1<int> px(2);
    //     vector1<int> py(2);

    //     for (int j = 0; j < (range2 - range1) - timeoffset1; j++)
    //     {
    //         for (int i = 0; i < NR - spaceoffset1; i++)
    //         {
    //             for (int k = 0; k < NC - spaceoffset1; k++)
    //             {
    //                 // cout << j << " ,(" << i << "," << k << ")" << endl;
    //                 bool gx = alldata[j](i, k);

    //                 bool gx2 = alldata[j + timeoffset1](periodicF(i + spaceoffset1, NR), k);

    //                 bool gx3 = alldata[j + timeoffset1](i, periodicF(k + spaceoffset1, NC));

    //                 pxy((int)gx, (int)gx2)++;

    //                 pxy((int)gx, (int)gx3)++;

    //                 px((int)(gx))++;

    //                 py((int)(gx2))++;
    //                 py((int)(gx3))++;
    //             }
    //         }
    // }

    // long int totp = pxy(0,0)+pxy(1,0)+pxy(0,1)+pxy(1,1);
    // long int totp2 = px[0]+px[1];
    // long int totp3 = py[0]+py[1];

    // double mi = 0.0;

   

    // for(int i1 = 0 ; i1 < 2 ; i1++) {
    //     for(int j1  = 0 ; j1 < 2 ; j1++) {
    //         double c1 = (double)pxy(i1,j1)/(double)(totp);

    //         double c2 = (double)px[i1] / (double)(totp2);

    //         double c3 = (double)py[j1] / (double)(totp3);
            
    //         // cout << c1 << " " << c2 << " " << c3 << endl;
    //         if(c1<1E-10) {

    //         mi+=0.;
    //         }
    //         else{

    //             mi += c1 * log2(c1 / (c2 * c3));
    //         }
    //     }
    // }
    
    // // cout << timeoffset1 << " " << spaceoffset1 << endl;
    // // cout << pxy << endl;
    // // cout << px << endl;
    // // cout << py << endl;
    // // cout << mic << endl;
    // // cout << mi << endl;

    // // cout << endl;

    // mia(timeoffset1,spaceoffset1) = mi;

    // }
    // }

    //cout << mia << endl;

    

    
}