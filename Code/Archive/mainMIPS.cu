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
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/unique.h>
#include <thrust/device_delete.h>
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
#include <curand.h>
#include <curand_kernel.h>


#include "DataStructures/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"


// #include "BrownianGel.cpp"
// #include "BrownianGel2.cpp"
// #include "LangevinGel.cpp"
// #include "LangevinGelFixed.cpp"

// #include "NCGasR.h"
#include "DensityFunctionalTheory/mipsGPU.cu"



using namespace std;

int main(int argc, char** argv) {
srand (time(NULL));


    int Nr = 200;
    int Nq = 40;

	float3 *d_pos;

	int size2 = Nr*Nr*Nq*sizeof(float3);


    cudaMalloc((void**)&d_pos,size2);

    float *d_xpos;
    float *d_qpos;

    int size3 = Nr*sizeof(float);
    int size4 = Nq*sizeof(float);

    cudaMalloc((void**)&d_xpos,size3);
    cudaMalloc((void**)&d_qpos,size4);
    
    int n = 1000;
    float l = 40.0;

    float volume_element = l*l*2.*pi/(float(Nr)*float(Nr)*float(Nq));

    assign_array<<<Nr,1>>>(d_xpos, l, Nr);
    assign_array<<<Nr,1>>>(d_qpos, 2*pi, Nq);

    assign_pos<<<Nr*Nr,Nq>>>(d_pos, d_xpos, d_qpos, Nr,Nq);

    potential_param p;
    p.l = l;
    p.epsilon = 1.0;
    p.sigma = 1.0;
    p.v0 = 1.0;

    float *d_den;
    cudaMalloc((void**)&d_den,Nr*Nr*Nq*sizeof(float));

    float val = (float)n/(l*l*2*pi);
    assign_den<<<Nr*Nr,Nq>>>(d_den,val,Nr*Nr*Nq);



    float *d_int;
    cudaMalloc((void**)&d_int,Nr*Nr*Nq*sizeof(float));
    for(int i = 0 ; i < Nr*Nr*Nq ; i++) {
    convo<<<Nr*Nr,Nq>>>(d_den, 
    d_pos,
    i,
    d_int,
    p,
    volume_element, Nr*Nr*Nq); 

    thrust::device_ptr<float> t_int(d_int);
   float tot = thrust::reduce(t_int,t_int+Nr*Nr*Nq);

    cout << i << " " << tot <<endl;
    }

return 0;

}