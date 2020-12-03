typedef struct {
    float l;
    float epsilon;
    float sigma;
    float v0;
} potential_param;

__device__ float potential(float3 r1, float3 r2, potential_param p) {
    float dx = r1.x - r2.x;
    float dy = r1.y - r2.y;
    float atemp = dx > 0 ? 1. : -1.;
    if(dx*dx > 0.25*p.l*p.l) {
        dx = dx - atemp*p.l;
    }
    float btemp = dy > 0 ? 1. : -1.;
    if(dy*dy > 0.25*p.l*p.l) {
        dy = dy - btemp*p.l;
    }

    float rSQR = SQR(dx) + SQR(dy);

    //float fac1 = (sigma / sqrt(rSQR));
    float fac2 = SQR(p.sigma)/rSQR;
    //else expf = 0.0;

    float st;
    if (rSQR > 0.69 * SQR(p.sigma))
    { //cut off
        float f1 = ((4. * p.epsilon * p.v0)) * (4 * fac2 * fac2 - 2 * fac2); //long range potential
        st = f1 * (dx * cos(r1.z) - dx * cos(r2.z) + dy * sin(r1.z) - dy * sin(r2.z));
    }
    else
    {
        st = 0.0;
    }

    return st;
}


// __global__ void calc_ints(float *den, float *res, int tot) {
// 	//particle positions in dev_unc_pos
// 	//dev cell list is the list of all particles with their cell
// 	int global_index = threadIdx.x + blockIdx.x * blockDim.x;

// 	if(global_index < tot) {
// 		dev_cell_list[global_index] = nbox*floorf(nbox*(dev_unc_pos[global_index].x)/lcell) +  floorf(nbox*(dev_unc_pos[global_index].y)/lcell);
// 		//dev_cell_list[global_index].y =
//     }
    
// }

template <class Q>
void print_device_array(Q *array, int n) {
Q *temparray = new Q [n];
cudaMemcpy(temparray,array,n*sizeof(Q),cudaMemcpyDeviceToHost);
for(int i = 0 ; i < n ; i++) {
	cout << temparray[i] <<  ",";
}
cout << endl;
delete temparray;

}


__global__ void assign_array(float *my_array, float len, int tot) {
   int global_index = threadIdx.x + blockIdx.x * blockDim.x;
    if(global_index < tot) {
         float dl = len/(float)(tot);
         my_array[global_index] = global_index*dl;
     }
}

__global__  void assign_pos(float3 *pos, float *xpos, float *qpos, int Nr, int Nq) {
    int global_index = threadIdx.x + blockIdx.x * blockDim.x;
    if(global_index < Nr*Nr*Nq) {
        float3 a;
        int q = global_index/SQR(Nr);
        int temp = global_index - q*SQR(Nr);
        int i = temp/Nr;
        int j = temp % Nr;

        a.x = xpos[i];
        a.y = xpos[j];
        a.z = qpos[q];
        pos[global_index] = a;
    }
}

__global__ void assign_den(float* den, float val, int n) {
    int global_index = threadIdx.x + blockIdx.x * blockDim.x;
    if(global_index < n) {
        den[global_index] = val;
    }
}

__global__ void convo(float *den, 
    float3 *points_in_space,
    int con,
    float *res,
    potential_param p,
    float volume_element, int tot) 
    {
    int global_index = threadIdx.x + blockIdx.x * blockDim.x;

    if(global_index < tot) {
        //int q = global_index % Nr

        res[global_index] =  volume_element* den[global_index] * potential(points_in_space[con],points_in_space[global_index],p);
    }
}


__global__ void calcint(float *den, float3 *points_in_space, float *res, potential_param p, float volume_element, int tot) {
    // int global_index = threadIdx.x + blockIdx.x * blockDim.x;
    // if(global_index < tot) {
    //     float *newres = new float [tot];

    //     convo(den,points_in_space,global_index,newres,p,volume_element,tot);

    //     thrust::device_ptr<float> t_int(newres);
    //     res[global_index] = thrust::reduce(t_int,t_int+tot);
        

        
    // }

}