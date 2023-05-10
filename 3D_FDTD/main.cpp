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
#define FREQ 10e9 // used for scaling simulation parameters (DT, DX, etc.)
#define PROP_FREQ FREQ // used for source wave

#define FRAC 20. // fraction of wavelength for spatial step (spatial resolution of grid)
#define CNO (0.7071) // Courant Number - stability factor


/* calculate useful values */
#define LAMBDA (C/FREQ) // wavelength
#define OMEGA (2*M_PI*FREQ)
#define PROP_OMEGA (2*M_PI*PROP_FREQ)
//#define BETA (OMEGA/C) // phase constant

/* define bounds of our simulation */
#define DX (LAMBDA/FRAC)
#define DY DX
#define DZ DX
#define DT (DX/C*CNO)
//#define X_MAX 100
//#define Y_MAX X_MAX
//#define Z_MAX Y_MAX
//#define N (int)ceil(X_MAX/DX)

#define T_MAX 200 // how long the simulation will run
//#define N (int)(FRAC * 3 / 2) // spatial size of grid
#define N 101
#define NX N
#define NY NX
#define NZ NX

//const int T_MAX, N, NX, NY;

// Create a struct template for E,H fields (to contain x,y,z,cond components)
struct field {
    double x[NX][NY][NZ]; // x component (x,y,z)
    double y[NX][NY][NZ]; // y component (x,y,z)
    double z[NX][NY][NZ]; // z component (x,y,z)

};

// Create a struct template for J,M fields (to contain x,y,z components)
struct current {
    double x[NX]; // x component
    double y[NX]; // y component
    double z[NX]; // z component

};

// Parameters of the simulation
struct parameters {
    double eps[NX][NY][NZ]; // relative permittivity of medium
    double mu[NX][NY][NZ]; // relative permeability of medium
    double ec[NX][NY][NZ]; // electrical conductivity of medium
    double mc[NX][NY][NZ]; // magnetic conductivity of medium
    int src_x; // x coordinate of source
    int src_y; // y coordinate of source
    int src_z; // z coordinate of source
    std::string test_name; // name of the output .txt file
};

/* function headers */
//int compute(field *E, field *H, current *J, current *M, parameters *params);
int compute(field *E, field *H, parameters *params);

int main() {
    std::cout << "Initializing Simulation..." << std::endl;
    std::cout << "NX: " << NX << " | NY: " << NY << "NZ: " << NZ << " | TMAX: " << T_MAX;
    std::cout << " | DX: " << DX << " | DY: " << DY << " | DZ: " << DZ << " | DT: " << DT << std::endl;

//  Initialize fields
    field E = {{0.,0., 0.},{0.,0., 0. },{0.,0., 0. }};
    field H = {{0.,0., 0.},{0.,0., 0. },{0.,0., 0. }};
//    current J = {{0.},{0.},{0.}};
//    current M = {{0.},{0.},{0.}};
    current *J = NULL;
    current *M = NULL;
    parameters  params = {{0., 0., 0.},{0.,0., 0.},
                          {0., 0., 0.},{0., 0., 0.},
                          0, 0, 0, ""};

    // adjust name of simulation output file
    params.test_name = "3D_test";

    params.src_x = round(NX/2);
    params.src_y = round(NY/2);
    params.src_z = round(NZ/2);


    int i, j, k;
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            for (k = 0; k < NZ; k++) {
                params.eps[i][j][k] = 1;
                params.mu[i][j][k] = 1;
                params.ec[i][j][k] = 0;
                params.mc[i][j][k] = 0;
//            if (i > NX / 2 && i < (NX/2 + (int)(LAMBDA/DX))) {
//                params.eps[i][j] = 1;
//            }
//                // rectangular waveguide - define conductive boundaries
//                int x1 = params.src_x - +round(LAMBDA / 4 / DY); // starting x location
////            x1 = params.src_x - 2;
////            int x2 = NX - 1; // ending x location
////            int wav_wid = 1.4 * LAMBDA / 2 / DY; // make width >lambda/2 to avoid cutoff
//                int wav_wid = 22.86e-3 / DY; // width of WR90 waveguide at 10GHz frequency
//                int y1 = params.src_y - floor(wav_wid / 2); // min y value (bottom)
//                int y2 = params.src_y + ceil(wav_wid / 2); // max y value (top)
//                if ((i > floor(x1 - 1) && i < floor(x1 + 1) && j > floor(y1 - 1) && j < floor(y2 + 1))
//                    //                    || (i > floor(x2 - 1) && i < floor(x2 + 1) && j > floor(y1 - 1) && j < floor(y2 + 1))
//                    || (j > floor(y1 - 1) && j < floor(y1 + 1) && i > floor(x1 - 1))
//                    || (j > floor(y2 - 1) && j < floor(y2 + 1) && i > floor(x1 - 1))
//                        ) {
//                    params.ec[i][j][k] = 1e9;
//                }
//                // rectangular waveguide - set dielectric thickness within metal
//                int x3 = round(NX / 2);
//                int sample_thickness = 1.1 * round(LAMBDA / 2 / DY); // multiply by factor of 1.1
//
//                if (i > x3 && i < x3 + sample_thickness && j > y1 - 1 && j < y2 + 1)
////                params.eps[i][j] = 10;
//                    ;
            }
        }
    }

    // run the simulation and store in file
//    compute(&E, &H, &J, &M, &params);
    int status = 0;

    status = compute(&E, &H, &params);

    if (status != 0)
        std::cout << "Warning: potential error in simulation!" << std::endl;

    return status;
}

//int compute(field *E, field *H, current *J, current *M, parameters *params) {
int compute(field *E, field *H, parameters *params) {
    // time-stepping loop
    std::cout << "Starting Simulation..." << std::endl;
    int t;
    std::string output_file = "/Users/williampoland/FDTD/2D/raw/" + params->test_name + ".txt";
    std::ofstream myfile (output_file);
    myfile << "FREQ T_MAX NX NY DT DX DY" << std::endl;
    myfile << PROP_FREQ << " " << T_MAX << " " << NX << " " << NY << " " << DT << " " << DX << " " << DY << std::endl;
//    myfile << "(t, x, y, z, Ez, Hx, Hy, eps_r, mu_r, e_cond, m_cond)" << std::endl;
    myfile << "(t, x, y, Ez, eps_r, ec)" << std::endl;

    // initialize variables for boundary conditions
    // these will store the past values of for edges of grid
    // 1st index - displacement from edge of boundary
    // 2nd index - previous time displacement
    // 3rd index - location along boundary
    const int num_disp = 3, time_disp = 2;
    double Ez_old_up[num_disp][time_disp][NX] = {0.,0.}, Ez_old_down[num_disp][time_disp][NX] = {0.,0.},
        Ez_old_left[num_disp][time_disp][NY] = {0.,0.}, Ez_old_right[num_disp][time_disp][NY] ={0.,0.};

    for (t = 0; t < T_MAX; t++) {
        int x, y;
        // step 1 - [TMz] Hx, Hy, Ez

        // Update Hx
        double cf_Hx, C_Hx_H, C_Hx_E; // coefficients for Hx
        for (x = 0; x < NX ; x++) {
            for (y = 0; y < NY - 1; y++) {
                cf_Hx = params->mc[x][y] * DT / 2 / MU_0 / params->mu[x][y];
                C_Hx_H = (1 - cf_Hx) / (1 + cf_Hx);
                C_Hx_E =  DT / MU_0 / params->mu[x][y] / DX / (1 + cf_Hx);

                H->x[x][y] = C_Hx_H * H->x[x][y] - C_Hx_E * (E->z[x][y+1] - E->z[x][y]);
            }
        }

        // Update Hy
        double cf_Hy, C_Hy_H, C_Hy_E; // coefficients for Hy
        for (x = 0; x < NX - 1 ; x++) {
            for (y = 0; y < NY; y++) {
                cf_Hy = params->mc[x][y] * DT / 2 / MU_0 / params->mu[x][y];
                C_Hy_H = (1 - cf_Hy) / (1 + cf_Hy);
                C_Hy_E = DT / MU_0 / params->mu[x][y] / DY / (1 + cf_Hy);

                H->y[x][y] = C_Hy_H * H->y[x][y] + C_Hy_E * (E->z[x+1][y] - E->z[x][y]);
            }
        }

        // Update Ez
        double cf_Ez, C_Ez_E, C_Ez_H;
        for (x = 1; x < NX - 1; x++) {
            for (y = 1; y < NY - 1; y++) {
                cf_Ez = params->ec[x][y] * DT / 2 / EPS_0 / params->eps[x][y];
                C_Ez_E = (1 - cf_Ez) / (1 + cf_Ez);
                C_Ez_H = DT / EPS_0 / params->eps[x][y] / DY / (1 + cf_Ez);

                E->z[x][y] = C_Ez_E * E->z[x][y] + C_Ez_H * ((H->y[x][y] - H->y[x-1][y]) - (H->x[x][y] - H->x[x][y-1]));
            }
        }

        // step 2 - [TEz] Ex, Ey, Hz

        // Update Ex
        double cf_Ex, C_Ex_E, C_Ex_H;
        for (x = 0; x < NX - 1; x++) {
            for (y = 1; y < NY - 1; y++) {
                cf_Ex = params->ec[x][y] * DT / 2 / EPS_0 / params->eps[x][y];
                C_Ex_E = (1 - cf_Ex) / (1 + cf_Ex);
                C_Ex_H = DT / EPS_0 / params->eps[x][y] / DX / (1 + cf_Ex);

                E->x[x][y] = C_Ex_E * E->x[x][y] + C_Ex_H * (H->z[x][y] - H->z[x][y-1]);
            }
        }

        // Update Ey
        double cf_Ey, C_Ey_E, C_Ey_H;
        for (x = 1; x < NX - 1; x++) {
            for (y = 0; y < NY - 1; y++) {
                cf_Ey = params->ec[x][y] * DT / 2 / EPS_0 / params->eps[x][y];
                C_Ey_E = (1 - cf_Ey) / (1 + cf_Ey);
                C_Ey_H = DT / EPS_0 / params->eps[x][y] / DY / (1 + cf_Ex);

                E->y[x][y] = C_Ey_E * E->y[x][y] + C_Ey_H * (H->z[x][y] - H->z[x-1][y]);
            }
        }

        // Update Hz
        double cf_Hz, C_Hz_E, C_Hz_H;
        for (x = 0; x < NX - 1; x++) {
            for (y = 0; y < NY - 1; y++) {
                cf_Hz = params->mc[x][y] * DT / 2 / MU_0 / params->mu[x][y];
                C_Hz_H = (1 - cf_Hz) / (1 + cf_Hz);
                C_Hz_E = DT / MU_0 / params->mu[x][y] / DX / (1 + cf_Ex);

                H->z[x][y] = C_Hz_H * H->z[x][y] + C_Hz_E * ((E->x[x][y+1] - E->x[x][y]) - (E->y[x+1][y] - E->y[x][y]));
            }
        }

        /* create source node */
        // Note: source nodes cannot be on boundaries x=0,x=NX,y=0,y=NY because of the absorbing boundaries
        // Source nodes for grid edges should be placed one node from the edge of grid (i.e. x=1 or y=NY-1)
        int delay = 0;
        int width = N/16;
        // hard-wired source
//        E->z[1][1] = exp(-(t - 30.) * (t - 30.) / 100.);
//        E->z[NX/2][NY/2] = 1.5 * exp(-(t - delay) * (t - delay) / width);
//            E->z[NX/2][NY/2] = 1;
//        E->z[NX/2][NY/2] = sin(CNO*BETA*t);
        // source for waveguide sim
        E->z[params->src_x][params->src_y] = sin(PROP_OMEGA*t*DT);
//        E->z[params->src_x][params->src_y] = sin(t);
        // vertical plane wave
//        for (int i = 0; i < NX; i++) {
////            E->z[i][0] = sin(t);
//            E->z[i][1] = exp(-(t - delay) * (t - delay) / width);
//        }
        // horizontal plane wave
//        for (int j = 0; j < NY; j++) {
//            E->z[1][j] = sin(0.5 * t);
////           E->z[1][j] = 1.5 * exp(-(t - delay) * (t - delay) / width);
//        }
        // additive source
//        E->z[NX/2][NY/2] += 1.5*exp(-(t - delay) * (t - delay) / width);




        /* second order ABCs for Ex,Ey,Ez @ x = 0, x = NX, y = 0, y = NY */

        // calculate coefficients
//        cf_Ez = params->ec[0][0] * DT / 2 / EPS_0 / params->eps[0][0];
//        C_Ez_H = DT / EPS_0 / params->eps[0][0] / DY / (1 + cf_Ez);
//
//        cf_Hy = params->mc[0][0] * DT / 2 / MU_0 / params->mu[0][0];
//        C_Hy_E = DT / MU_0 / params->mu[0][0] / DY / (1 + cf_Hy);
//
//        double temp1_TM = sqrt(C_Ez_H * C_Hy_E);
//        double temp2_TM = 1. / temp1_TM + 2. + temp1_TM;
//        double cf0 = -(1. / temp1_TM - 2. + temp1_TM) / temp2_TM;
//        double cf1 = -2. * (temp1_TM - 1. / temp1_TM) / temp2_TM;
//        double cf2 = 4. * (temp1_TM + 1. / temp1_TM) / temp2_TM;
//
//
//        // ABCs @ x=0 & x=NX
//        for (y = 0; y < NY; y++) {
//            // calculate ABC @ x=0
////            E->y[0][y] = 0;
//
//            E->z[0][y] = cf0 * (E->z[2][y] + Ez_old_left[0][1][y])
//                + cf1 * (Ez_old_left[0][0][y] + Ez_old_left[2][0][y] - E->z[1][y] - Ez_old_left[1][1][y])
//                + cf2 * Ez_old_left[1][0][y]
//                - Ez_old_left[2][1][y];
//            // update old values @ x=0
//            for (x = 0; x < 3; x++) {
//                Ez_old_left[x][1][y] = Ez_old_left[x][0][y];
//                Ez_old_left[x][0][y] = E->z[x][y];
//            }
//
//            // calculate ABC @ x=NX
////            E->x[NX - 1][y] = 0;
//
//            E->z[NX - 1][y] = cf0 * (E->z[NX - 3][y] + Ez_old_right[0][1][y])
//                + cf1 * (Ez_old_right[0][0][y] + Ez_old_right[2][0][y] - E->z[NX - 2][y] - Ez_old_right[1][1][y])
//                + cf2 * Ez_old_right[1][0][y]
//                - Ez_old_right[2][1][y];
//            // update old values @ x=0
//            for (x = 0; x < 3; x++) {
//                Ez_old_right[x][1][y] = Ez_old_right[x][0][y];
//                Ez_old_right[x][0][y] = E->z[NX - 1 - x][y];
//            }
//        }
//        // ABCs @ y=0 & y=NY
//        for (x = 0; x <NX; x++) {
//            // calculate ABC @ y=0
////            E->x[x][0] = 0;
//
//            E->z[x][0] = cf0 * (E->z[x][2] + Ez_old_down[0][1][x])
//                + cf1 * (Ez_old_down[0][0][x] + Ez_old_down[2][0][x] - E->z[x][1] - Ez_old_down[1][1][x])
//                + cf2 * Ez_old_down[1][0][x]
//                - Ez_old_down[2][1][x];
//            // update old values @ y=0
//            for (y = 0; y < 3; y++) {
//                Ez_old_down[y][1][x] = Ez_old_down[y][0][x];
//                Ez_old_down[y][0][x] = E->z[x][y];
//            }
//
//            // calculate ABC @ y=NY
////            E->x[x][NY - 1] = 0;
//
//            E->z[x][NY - 1] = cf0 * (E->z[x][NY - 3] + Ez_old_up[0][1][x])
//                + cf1 * (Ez_old_up[0][0][x] + Ez_old_up[2][0][x] - E->z[x][NY - 2] - Ez_old_up[1][1][x])
//                + cf2 * Ez_old_up[1][0][x]
//                - Ez_old_up[2][1][x];
//            // update old values @ y=0
//            for (y = 0; y < 3; y++) {
//                Ez_old_up[y][1][x] = Ez_old_up[y][0][x];
//                Ez_old_up[y][0][x] = E->z[x][NY - 1 - y];
//            }
//        }





//        for (int)

        // now write relevant data to file
        if (myfile.is_open())
        {
            std::cout << "Writing to file...[t=" << t << "]" << std::endl;
            for (x = 0; x < NX; x++) {
                for (y = 0; y <NY; y++) {
                    myfile << t << " " << x << " " << y << " " << E->z[x][y] << " " << params->eps[x][y] << " " << params->ec[x][y] << std::endl;
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

    return 0;
}

