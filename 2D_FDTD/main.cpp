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

#define FRAC 20. // fraction of wavelength for spatial step (spatial resolution of grid)
#define CNO (0.7071) // Courant Number - stability factor


/* calculate useful values */
#define LAMBDA (C/FREQ) // wavelength
#define BETA (2*M_PI*FREQ/C) // phase constant

/* define bounds of our simulation */
#define DX (LAMBDA/FRAC)
#define DY DX
#define DT (DX/C*CNO)
//#define X_MAX 100
//#define Y_MAX X_MAX
//#define Z_MAX Y_MAX
//#define N (int)ceil(X_MAX/DX)

#define T_MAX 100 // how long the simulation will run
//#define N (int)(FRAC * 3 / 2) // spatial size of grid
#define N 100
#define NX N
#define NY NX


// Create a struct template for E,H fields (to contain x,y,z,cond components)
struct field {
    double x[NX][NY]; // x component (x,y)
    double y[NX][NY]; // y component (x,y)
    double z[NX][NY]; // z component (x,y)

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
    std::cout << "Initializing Simulation..." << std::endl;
    std::cout << "NX: " << NX << " | NY: " << NY << " | TMAX: " << T_MAX;
    std::cout << " | DX: " << DX << " | DY: " << DY << " | DT: " << DT << std::endl;

//  Initialize fields
    field E = {{0.,0.},{0.,0. },{0.,0. }};
    field H = {{0.,0.},{0.,0. },{0.,0. }};
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
    std::cout << "Starting Simulation..." << std::endl;
    int t;
    std::string test_name = "basic_sine1Hz";
    std::string output_file = "/Users/williampoland/FDTD/2D/raw/" + test_name + ".txt";
    std::ofstream myfile (output_file);
    myfile << "T_MAX; NX; NY" << std::endl;
    myfile << T_MAX << std::endl << NX << std::endl << NY << std::endl;
//    myfile << "(t, x, y, z, Ez, Hx, Hy, eps_r, mu_r, e_cond, m_cond)" << std::endl;
    myfile << "(t, x, y, Ez)" << std::endl;

    // initialize variables for boundary conditions
//    double Ez1_old = 0; // we will use this to store past values of E->z[1]
//    double EzNm2_old = 0; //E->z[N-2]

    for (t = 0; t < T_MAX; t++) {
        int x, y;
        // step 1 - [TMz] Hx, Hy, Ez

        // Update HX
        double C_Hx_H, C_Hx_E; // coefficients for Hx
        for (x = 0; x < NX ; x++) {
            for (y = 0; y < NY - 1; y++) {
                C_Hx_H = 1;
                C_Hx_E =  DT / MU_0 / DX;
                H->x[x][y] += -C_Hx_E * (E->z[x][y+1] - E->z[x][y]);
//                H->x[x][y] = (1. - C_Hx) / (1. + C_Hx) * H->x[x][y]
//                        - DT / (1. + C_Hx) / DY / MU_0 / params->mu[x] * (E->z[x][y+1] - E->z[x][y]);
            }
        }

        // Update HY
        double C_Hy_H, C_Hy_E; // coefficients for Hy
        for (x = 0; x < NX - 1 ; x++) {
            for (y = 0; y < NY; y++) {
                C_Hy_H = 1;
                C_Hy_E = DT / MU_0 / DY;
                H->y[x][y] += C_Hy_E * (E->z[x+1][y] - E->z[x][y]);
//                H->y[x][y] = (1. - C_Hy) / (1. + C_Hy) * H->y[x][y]
//                        + DT / (1. + C_Hy) / DX / MU_0 / params->mu[x] * (E->z[x][y+1] - E->z[x][y]);
            }
        }

        // Update EZ
        double C_Ez_E, C_Ez_H;
        for (x = 1; x < NX - 1; x++) {
            for (y = 1; y < NY - 1; y++) {
                C_Ez_E = 1;
                C_Ez_H = DT / EPS_0 / DX;
                E->z[x][y] += C_Ez_H * ((H->y[x][y] - H->y[x-1][y]) - (H->x[x][y] - H->x[x][y-1]));
//                E->z[x][y] = (1. - C_Ez) / (1. + C_Ez) * E->z[x][y]
//                        + DT / (1. + C_Ez) / params->eps[x] * (1. / DX / EPS_0 / params->eps[x] * (H->y[x][y] - H->y[x-1][y])
//                        - 1. / DY / EPS_0 / params->eps[x] * (H->x[x][y] - H->x[x][y-1]));
            }
        }

        /* create source node */
        int delay = 0;
        int width = N/4;
        // hard-wired source
//        E->z[0][0] = exp(-(t - 30.) * (t - 30.) / 100.);
//        E->z[NX/2][NY/2] = exp(-(t - delay) * (t - delay) / width);
//            E->z[NX/2][NY/2] = 1;
//        E->z[NX/2][NY/2] = sin(CNO*BETA*t);
        E->z[NX/2][NY/2] = sin(1*t);
        // additive source
//        E->z[NX/2][NY/2] += exp(-(t - delay) * (t - delay) / width);





        /* boundary conditions*/
        // 1st order ABC at x = 0
//        double coef = C * DT / DX / sqrt(params->eps[0]*params->mu[0]);
//        E->z[0] = Ez1_old + (coef - 1) / (coef + 1) * (E->z[1]- E->z[0]); // update boundaries
//        E->z[N-1] = EzNm2_old + (coef - 1) / (coef + 1) * (E->z[N-2]- E->z[N-1]);
//        Ez1_old = E->z[1]; // update past values
//        EzNm2_old = E->z[N-2];

        // now write relevant data to file
        if (myfile.is_open())
        {
            std::cout << "Writing to file...[t=" << t << "]" << std::endl;
            for (x = 0; x < NX; x++) {
                for (y = 0; y <NY; y++) {
                    myfile << t << " " << x << " " << y << " " << E->z[x][y] << std::endl;
//                    myfile << t << " " << x << " " << y << " "
//                        << E->z[x][y] << " " << H->x[x][y] << " " << H->y[x][y] << " "
//                        << params->eps[x] << " " << params->mu[x] << " " << params->ec[x] << " " << params->mc[x]
//                        << std::endl;
                }
            }
        }
        else std::cout << "Unable to open file" << std::endl;
    }
    myfile.close();
    std::cout << "Finished Writing to File: " << output_file << std::endl;
    std::cout << "Total entries: " << NX * NY * T_MAX << std::endl;
    std::cout << "Simulation Complete." << std::endl;
}

