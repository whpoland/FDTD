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

#define T_MAX 500
#define N 400 // the number of elements we will track at a time


// Create a struct template for E,H fields (to contain x,y,z,cond components)
struct field {
    double x[N]; // x component
    double y[N]; // y component
    double z[N]; // z component

};

// Create a struct template for J,M fields (to contain x,y,z components)
struct current {
    double x[N]; // x component
    double y[N]; // y component
    double z[N]; // z component

};

// Parameters of the medium
struct parameters {
    double eps[N]; // relative permittivity
    double mu[N]; // relative permeability
    double ec[N]; // electrical conductivity
    double mc[N]; // magnetic conductivity
};

/* function headers */
void compute(field *E, field *H, current *J, current *M, parameters *params);

int main() {
//  Initialize fields
    field E = {{0.},{0.},{0.}};
    field H = {{0.},{0.},{0.}};
    current J = {{0.},{0.},{0.}};
    current M = {{0.},{0.},{0.}};
    parameters  params = {{0.},{0.},{0.},{0.}};

    for (int i = 0; i < N; i++) {
        params.eps[i] = 1;
        params.mu[i] = 1;
        params.ec[i] = 0;
        params.mc[i] = 0;
        if (i < N/2) {
//            params.eps[i] = 1;
//            params.ec[i] = 20;
        }
        else {
//            params.eps[i] = 10;
        }
    }

    // run the simulation and store in file
    compute(&E, &H, &J, &M, &params);

    return 0;
}

void compute(field *E, field *H, current *J, current *M, parameters *params) {
    // time-stepping loop
    int t;
    std::string output_file = "/Users/williampoland/FDTD/1D/midsrc_ABC.txt";
    std::ofstream myfile (output_file);

    // initialize variables for boundary conditions
    double Ez1_old = 0; // we will use this to store past values of E->z[1]
    double EzNm2_old = 0; //E->z[N-2]

    for (t = 0; t < T_MAX; t++) {
        // step 1 - [TE] Hx, Hy, Ez
        int x;
        for (x = 0; x < N - 1; x++) {
//            H->y[x] += (E->z[x+1] - E->z[x]) * DT / MU_0 / params->mu[x] / DX;
            H->y[x] += DT / MU_0 / params->mu[x] * ((E->z[x+1] - E->z[x]) / DX + H->y[x] * params->mc[x] + M->y[x]);
            //never reaches Hy[N-1], so this node stays 0 -- PMC
        }

        for (x = 1; x < N; x++) {
//            E->z[x] += (H->y[x] - H->y[x-1]) * DT / EPS_0 / params->eps[x] / DX;
            E->z[x] += DT / EPS_0 / params->eps[x] * ((H->y[x] - H->y[x-1]) / DX - E->z[x] * params->ec[x] - J->z[x]);
            //never reaches Ez[0], so this node stays 0 -- PEC
        }
        // step 2 - [TM] Ex, Ey, Hz


        /* hardwire a source node */
        // this makes Ez[0] follow a  wave shape, forcing the simulation to create a wave
//        E->z[0] = exp(-(t - 30.) * (t - 30.) / 100.);

        /* additive source node */
        // middle of grid
        E->z[N/2] += exp(-(t - 30.) * (t - 30.) / 100.);

        // end of grid
//        E->z[N - 1] += exp(-(t - 30.) * (t - 30.) / 100.);


        /* boundary conditions*/
        // 1st order ABC at x = 0
        double coef = C * DT / DX / sqrt(params->eps[0]*params->mu[0]);
        E->z[0] = Ez1_old + (coef - 1) / (coef + 1) * (E->z[1]- E->z[0]); // update boundaries
        E->z[N-1] = EzNm2_old + (coef - 1) / (coef + 1) * (E->z[N-2]- E->z[N-1]);
        Ez1_old = E->z[1]; // update past values
        EzNm2_old = E->z[N-2];

        // 2nd order ABC at x = 0
//        double coef1 = C * DT / DX;
//        double coef2 = coef1 / sqrt(params->eps[0]*params->mu[0]);
//        E->z[0] = -1 / (1/coef2 + 2 + coef2) * ((1/coef2 - 2 + coef2)*(E->z[2]+E->z[0])+2*(coef2 - 1/coef2)*() - 4*(1/coef2 +coef2)*E->z[1])- E->z[2];

        // now write relevant data to file
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

