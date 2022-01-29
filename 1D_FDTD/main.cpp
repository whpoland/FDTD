#include <iostream>
#include <fstream>
#include <cmath>
//#include "matplotlib-cpp-master/matplotlibcpp.h"

//namespace plt = matplotlibcpp;

/* define mathematical constants */
#define C 3e8
#define EPS_0 8.85e-12
#define MU_0 (M_PI*4e-7)

/* CUSTOMIZABLE VALUES */
#define FREQ 300e6

/* define bounds of our simulation */
#define LAMBDA C/FREQ

#define DX (LAMBDA/10.)
#define DT (DX/C)
//#define X_MAX 100
//#define Y_MAX X_MAX
//#define Z_MAX Y_MAX
//#define N (int)ceil(X_MAX/DX)

#define T_MAX 800
#define N 400 // the number of elements we will track at a time


// Create a struct template for E,H fields (to contain x,y,z,cond components)
struct field {
    double x[N]; // x component
    double y[N]; // y component
    double z[N]; // z component
    double c[N]; // conductivity
};

// Create a struct template for J,M fields (to contain x,y,z components)
struct current {
    double x[N]; // x component
    double y[N]; // y component
    double z[N]; // z component
};

struct parameters {
    double eps[N]; // relative permittivity
    double mu[N]; // relative permeability
    double ec[N]; // electrical conductivity
    double mc[N]; // magnetic conductivity
};

/* function headers */
double numDeriv(double* arr, int i, int delta);
double updateEx(int t, int delta, field* E, field *H, current* J, current* M);
double updateEy(int t, int delta, field* E, field *H, current* J, current* M);
double updateEz(int t, int delta, field* E, field *H, current* J, current* M);
double updateHx(int t, int delta, field* E, field *H, current* J, current* M);
double updateHy(int t, int delta, field* E, field *H, current* J, current* M);
double updateHz(int t, int delta, field* E, field *H, current* J, current* M);
void compute(field *E, field *H, current *J, current *M, parameters *params);

int main() {
//  Initialize fields
    field E = {{0.},{0.},{0.},{0.}};
    field H = {{0.},{0.},{0.},{0.}};
    current J = {{0.},{0.},{0.}};
    current M = {{0.},{0.},{0.}};
    parameters  params = {{0.},{0.},{0.},{0.}};

    for (int i = 0; i < N; i++) {
        params.eps[i] = 1;
        params.mu[i] = 1;
        params.ec[i] = 0;
        params.mc[i] = 0;
        if (i < N/2) {
            params.eps[i] = 10;
        }
    }

//    params.ec[100] = 100;

//    std::cout << std::endl;
//    std::cout << "EPS0: " << EPS_0 << " MU0: " << MU_0 << " DT: " << DT << " DX: " << DX << std::endl;
//    std::cout << "DT/DX: " << DT / DX << std::endl;
//    std::cout << "LAMBDA: " << LAMBDA << std::endl;
//    std::cout << "C: " << C << std::endl;
//    std::cout << "DT/MU/DX: " << DT / MU_0 / DX << std::endl;
//    std::cout << "DT/EPS/DX: " << DT / EPS_0 / DX << std::endl;


    // run the simulation and store in file
    compute(&E, &H, &J, &M, &params);

    return 0;
}

double numDeriv(double* arr, int i, int delta) {
    return (arr[i + delta/2] + arr[i - delta/2])/delta;
}

double updateEx(int t, int delta, field* E, field *H, current* J, current* M) {

}

double updateHx(int t, int delta, field* E, field *H, current* J, current* M) {
//    H->x[t - delta/2] + delta/MU*(numDeriv(E->z, t, delta) + M->x[t] + H->c[t]*H->x[t])
}

void compute(field *E, field *H, current *J, current *M, parameters *params) {
    // time-stepping loop
    int t;
    std::string output_file = "/Users/williampoland/FDTD Output/FDTDdie.txt";
//    std::string output_file = "/Users/williampoland/FDTD Output/FDTD3.txt";
//    std::string output_file = "/Users/williampoland/FDTD Output/FDTDadd.txt";
    std::ofstream myfile (output_file);

    for (t = 0; t < T_MAX; t++) {
        // step 1 - [TE] Hx, Hy, Ez
        int x;
        for (x = 0; x < N - 1; x++) {
//            H->y[x] += (E->z[x + 1] - E->z[x]) * DT / MU_0 / params->mu[x] / DX;
//            H->y[x] += (E->z[x+1] - E->z[x]) * DT / MU_0 / params->mu[x] / DX;
            H->y[x] += DT / MU_0 / params->mu[x] * ((E->z[x+1] - E->z[x]) / DX + H->y[x] * params->mc[x] + M->y[x]);
            //never reaches Hy[N-1], so this node stays 0 -- PMC
        }
        // step 2 - [TM] Ex, Ey, Hz
        for (x = 1; x < N; x++) {
//            E->z[x] += (H->y[x] - H->y[x - 1]) * DT / EPS_0 / DX;
//            E->z[x] += (H->y[x] - H->y[x-1]) * DT / EPS_0 / params->eps[x] / DX;
            E->z[x] += DT / EPS_0 / params->eps[x] * ((H->y[x] - H->y[x-1]) / DX - E->z[x] * params->ec[x] - J->z[x]);
        }
        /* hardwire a source node */
        // this makes Ez[0] follow a  wave shape, forcing the simulation to create a wave
//        E->z[0] = exp(-(t - 30.) * (t - 30.) / 100.);

        /* additive source node */
        E->z[N - 1] += exp(-(t - 30.) * (t - 30.) / 100.);

        if (myfile.is_open())
        {
            for (x = 0; x < N; x++)
                myfile << E->z[x] << " ";
        }
        else std::cout << "Unable to open file" << std::endl;

        myfile << "\n";
//        printf("%d\n",E->z[50]);
    }
    myfile.close();
    std::cout << std::endl << "Finished Writing to File: " << output_file << std::endl;
}

