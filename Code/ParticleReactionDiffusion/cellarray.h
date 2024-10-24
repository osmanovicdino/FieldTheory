#ifndef CELLARRAY_H
#define CELLARRAY_H

#include "prd.h"


struct pp {
    int i;
    int j;
    double time;
};

struct self_interaction
{
    virtual self_interaction *clone() const = 0; // Virtual constructor (clone method)
    virtual double operator()(int N)=0; 
};
struct attractive_self_interaction : self_interaction {
    double eps;
    double Neq; //this is the optimal number for energy in N
    double operator()(int N) {
        return -2*eps*SQR(N)+(eps/SQR(Neq))* SQR(N*N);
    }
    attractive_self_interaction *clone() const
    {
        return new attractive_self_interaction(*this);
    }
};
struct no_self_attraction : self_interaction {
    double operator()(int N) {return 0.;}
    no_self_attraction *clone() const
    {
        return new no_self_attraction(*this);
    }
};

struct surface_interaction {
    virtual double operator()(int N1,int N2)=0;
    virtual surface_interaction *clone() const = 0; // Virtual constructor (clone method)
};

struct linear_surface_interaction : surface_interaction {
    double eps;
    double operator()(int N1, int N2) {
        return eps*N1*N2;
    }
    linear_surface_interaction *clone() const {
        return new linear_surface_interaction(*this);
    } // Virtual constructor (clone method)
};
struct no_surface_interaction : surface_interaction
{
    double operator()(int N1, int N2) {return 0.0;}
    no_surface_interaction *clone() const {
        return new no_surface_interaction(*this);
    } // Virtual constructor (clone method)
};


struct diffusionstore2 {
    diffusionstore a;
    int posi;
    int posj;
    bool operator<(const diffusionstore2 &other)
    {
        return a < other.a;
    }
};

struct cell_array {

int Ntypes;
vector1<self_interaction*> si;
matrix<surface_interaction *> cross_si;
matrix<surface_interaction *> cnni;

vector<diffusionstore2> next_diffusion_events;



double current_time;
int Nx;

matrix<cell> reactors; // each matrix is a cell, the cells are taken to be a well mixed chemical reactor

// pp next_diffusion_event;
// pp next_reaction_event;


cell_array(int, int , cell&);

void initialize_next_diffusion_events();

void set_si_interaction(const self_interaction &pot, int i) { si[i] = pot.clone(); }
void set_crosssi_interaction(const surface_interaction &pot, int i, int j) { cross_si(i,j) = pot.clone(); };
void set_cnni_interaction(const surface_interaction &pot, int i, int j) { cnni(i, j) = pot.clone(); }

double getenergy(int,int,int,int,int);

void find_next_diffusion_event();
void find_next_reaction_event();

void collision_check();

void do_a_diffusion_event(vector1<int>&);
void diffusionmove(int,int,int,int,int);
void diffusionstay(int,int,int);


void iterate();

void reset_all_times();
};


#include "cellarray.cpp"

#endif