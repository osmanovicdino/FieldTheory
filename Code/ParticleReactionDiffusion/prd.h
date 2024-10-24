#ifndef PRD_H
#define PRD_H

#include "../DataStructures/vector1.h"
#include "../DataStructures/matrix2.h"



struct Diffusion_Rule {

};

struct Subdiffusion {
bool normal;
double alpha;
double tau;
double p;

double generate_alarm() {
    if(normal) {
        double r = (double)rand() / (double)(RAND_MAX);
        return log(1. / (1. - r)) /tau;
    }
    else{
    double r = (double)rand() / (double)(RAND_MAX);
    return -(((-1 + pow(1 - r, 1. /alpha)) *tau) / pow(1 - r, 1. /alpha));
    }
}

double generate_initial() { //the initial condition is somewhat different from generating a new alarm time (see mathematica workbook)
    if(normal) {
        double r = (double)rand() / (double)(RAND_MAX);
        return log(1. / (1. - r)) / tau;
    }
    else{
    double r = (double)rand() / (double)(RAND_MAX);

    return -tau + pow(-(((p + tau - pow(tau, 1 - alpha) * pow(p + tau, alpha)) * (-r - (pow(tau, 1 - alpha) * pow(p + tau, alpha)) / (p + tau - pow(tau, 1 - alpha) * pow(p + tau, alpha)))) /
                 pow(p + tau, alpha)),
               1 / (1 - alpha));
    }
}

};

// struct Ndiffusion
// {
//     double lambda;
//     double generate_alarm()
//     {
//         double r = (double)rand() / (double)(RAND_MAX);
//         return log(1./(1.-r))/lambda;
//     }

//     double generate_initial()
//     { // the initial condition is somewhat different from generating a new alarm time (see mathematica workbook)

//         double r = (double)rand() / (double)(RAND_MAX);
//         return log(1. / (1. - r)) / lambda;
//     }
// };

struct chemical_reaction {
int orderin;
int orderout;
double rate;
int i1,i2;
int o1,o2;

};


struct diffusionstore {
    int nt;
    //int wh;
    double time;
    bool operator<(const diffusionstore &other) {
        return time < other.time;
    }

};
ostream& operator<<(ostream &s, const diffusionstore &a) {
    s << a.nt << "," <<  a.time;
    return s;
}
	

void insert_new(vector<diffusionstore> &prev, diffusionstore a)
{
    auto insertion_point = std::lower_bound(prev.begin(), prev.end(), a);
    prev.insert(insertion_point, a);
}

struct cell {
    
    int Ntypes; //the number of types of particles
    int Nreactions;
    Subdiffusion *sd; //the diffusion properties of each of the particles
    chemical_reaction *cr;

    //vector1<vector<double>> diffusion_times;
    vector<diffusionstore> sorted_times;
    vector1<int> typenumber; //keep track of the type number

    double reaction_clock;
    int which_reaction;
    double diffusion_clock;
    int which_diffusion_type;
    int which_diffusion;
    
    void initialize(vector1<int>&);
    void set_chemical_reaction(chemical_reaction, int);
    void set_subdiffusion(Subdiffusion, int);

    void find_random_single(int, int &);
    void find_random_pair(int,int,int&,int&);

    void generate_reaction_time(double);
    void do_reaction(double);

    void generate_diffusion_time();

    

    cell();
    cell(int a) : cell() {}
    cell(int Nt, int Nr);
    cell(const cell &a);
    cell &operator=(const cell &); // assignment operator
                                                           //  Gillespie Algorithm is a lot simpler
    ~cell() {
        delete sd;
        delete cr;
    }
    //when a reaction happens
};


#include "prd.cpp"

#endif