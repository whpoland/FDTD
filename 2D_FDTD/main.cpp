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
#define OMEGA (2*M_PI*FREQ)
#define BETA (OMEGA/C) // phase constant

/* define bounds of our simulation */
#define DX (LAMBDA/FRAC)
#define DY DX
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

//const int T_MAX, N, NX, NY;

// Create a struct template for E,H fields (to contain x,y,z,cond components)
struct field {
    double x[NX][NY]; // x component (x,y)
    double y[NX][NY]; // y component (x,y)
    double z[NX][NY]; // z component (x,y)

};

// Create a struct template for J,M fields (to contain x,y,z components)
struct current {
    double x[NX]; // x component
    double y[NX]; // y component
    double z[NX]; // z component

};

// Parameters of the medium
struct parameters {
    double eps[NX][NY]; // relative permittivity
    double mu[NX]; // relative permeability
    double ec[NX][NY]; // electrical conductivity
    double mc[NX]; // magnetic conductivity
    int src_x; // x coordinate of source
    int src_y; // y coordinate of source
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
    parameters  params = {{0., 0.},{0.},{0., 0.},{0.}, 0, 0};

    //    params.src_x = 2 + ceil(LAMBDA / 4 / DY);
    params.src_x = 3;
    params.src_y = round(NY/2);

    int i, j;
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            params.eps[i][j] = 1;
            params.mu[i] = 1;
            params.ec[i][j] = 0;
            params.mc[i] = 0;
//            if (i > NX / 2 && i < (NX/2 + (int)(LAMBDA/DX))) {
//                params.eps[i][j] = 4;
////                params.ec[i][j] = 1;
//            }
            // rectangular waveguide - define conductive boundaries
            int x1 = params.src_x - 1; // starting x location
//            int x2 = NX - 1; // ending x location
            int wav_wid = LAMBDA / 2 / DY;
            int y1 = params.src_y - floor(wav_wid/2); // min y value (bottom)
            int y2 = params.src_y + ceil(wav_wid/2); // max y value (top)
            if ((i > floor(x1 - 1) && i < floor(x1 + 1) && j > floor(y1 - 1) && j < floor(y2 + 1))
//                    || (i > floor(x2 - 1) && i < floor(x2 + 1) && j > floor(y1 - 1) && j < floor(y2 + 1))
                    || (j > floor(y1 - 1) && j < floor(y1 + 1) && i > floor(x1 - 1))
                    || (j > floor(y2 - 1) && j < floor(y2 + 1) && i > floor(x1 - 1))
                    ){
                params.ec[i][j] = 2;
            }
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
    std::string test_name = "waveguide3";
    std::string output_file = "/Users/williampoland/FDTD/2D/raw/" + test_name + ".txt";
    std::ofstream myfile (output_file);
    myfile << "FREQ T_MAX NX NY DT DX DY" << std::endl;
    myfile << FREQ << " " << T_MAX << " " << NX << " " << NY << " " << DT << " " << DX << " " << DY << std::endl;
//    myfile << "(t, x, y, z, Ez, Hx, Hy, eps_r, mu_r, e_cond, m_cond)" << std::endl;
    myfile << "(t, x, y, Ez, ec)" << std::endl;

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

        // Update HX
        double cf_Hx, C_Hx_H, C_Hx_E; // coefficients for Hx
        for (x = 0; x < NX ; x++) {
            for (y = 0; y < NY - 1; y++) {
                cf_Hx = params->mc[x] * DT / 2 / MU_0 / params->mu[x];
                C_Hx_H = (1 - cf_Hx) / (1 + cf_Hx);
                C_Hx_E =  DT / MU_0 / params->mu[x] / DX / (1 + cf_Hx);

                H->x[x][y] = C_Hx_H * H->x[x][y] - C_Hx_E * (E->z[x][y+1] - E->z[x][y]);
            }
        }

        // Update HY
        double cf_Hy, C_Hy_H, C_Hy_E; // coefficients for Hy
        for (x = 0; x < NX - 1 ; x++) {
            for (y = 0; y < NY; y++) {
                cf_Hy = params->mc[x] * DT / 2 / MU_0 / params->mu[x];
                C_Hy_H = (1 - cf_Hy) / (1 + cf_Hy);
                C_Hy_E = DT / MU_0 / params->mu[x] / DY / (1 + cf_Hy);

                H->y[x][y] = C_Hy_H * H->y[x][y] + C_Hy_E * (E->z[x+1][y] - E->z[x][y]);
            }
        }

        // Update EZ
        double cf_Ez, C_Ez_E, C_Ez_H;
        for (x = 1; x < NX - 1; x++) {
            for (y = 1; y < NY - 1; y++) {
                cf_Ez = params->ec[x][y] * DT / 2 / EPS_0 / params->eps[x][y];
                C_Ez_E = (1 - cf_Ez) / (1 + cf_Ez);
                C_Ez_H = DT / EPS_0 / params->eps[x][y] / DY / (1 + cf_Ez);

                E->z[x][y] = C_Ez_E * E->z[x][y] + C_Ez_H * ((H->y[x][y] - H->y[x-1][y]) - (H->x[x][y] - H->x[x][y-1]));
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
        E->z[params->src_x][params->src_y] = sin(OMEGA*t*DT);
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




        /* second order ABCs for x = 0, x = NX, y = 0, y = NY */
        // calculate coefficients
        cf_Ez = params->ec[0][0] * DT / 2 / EPS_0 / params->eps[0][0];
        C_Ez_H = DT / EPS_0 / params->eps[0][0] / DY / (1 + cf_Ez);

        cf_Hy = params->mc[0] * DT / 2 / MU_0 / params->mu[0];
        C_Hy_E = DT / MU_0 / params->mu[0] / DY / (1 + cf_Hy);

        double temp1 = sqrt(C_Ez_H * C_Hy_E);
        double temp2 = 1. / temp1 +  2. + temp1;
        double cf0 = -(1. / temp1 - 2. + temp1) / temp2;
        double cf1 = -2. * (temp1 - 1. / temp1) / temp2;
        double cf2 = 4. * (temp1 + 1. / temp1) / temp2;

        // ABCs @ x=0 & x=NX
        for (y = 0; y < NY; y++) {
            // calculate ABC @ x=0
            E->z[0][y] = cf0 * (E->z[2][y] + Ez_old_left[0][1][y])
                + cf1 * (Ez_old_left[0][0][y] + Ez_old_left[2][0][y] - E->z[1][y] - Ez_old_left[1][1][y])
                + cf2 * Ez_old_left[1][0][y]
                - Ez_old_left[2][1][y];
            // update old values @ x=0
            for (x = 0; x < 3; x++) {
                Ez_old_left[x][1][y] = Ez_old_left[x][0][y];
                Ez_old_left[x][0][y] = E->z[x][y];
            }

            // calculate ABC @ x=NX
            E->z[NX - 1][y] = cf0 * (E->z[NX - 3][y] + Ez_old_right[0][1][y])
                + cf1 * (Ez_old_right[0][0][y] + Ez_old_right[2][0][y] - E->z[NX - 2][y] - Ez_old_right[1][1][y])
                + cf2 * Ez_old_right[1][0][y]
                - Ez_old_right[2][1][y];
            // update old values @ x=0
            for (x = 0; x < 3; x++) {
                Ez_old_right[x][1][y] = Ez_old_right[x][0][y];
                Ez_old_right[x][0][y] = E->z[NX - 1 - x][y];
            }
        }
        // ABCs @ y=0 & y=NY
        for (x = 0; x <NX; x++) {
            // calculate ABC @ y=0
            E->z[x][0] = cf0 * (E->z[x][2] + Ez_old_down[0][1][x])
                + cf1 * (Ez_old_down[0][0][x] + Ez_old_down[2][0][x] - E->z[x][1] - Ez_old_down[1][1][x])
                + cf2 * Ez_old_down[1][0][x]
                - Ez_old_down[2][1][x];
            // update old values @ y=0
            for (y = 0; y < 3; y++) {
                Ez_old_down[y][1][x] = Ez_old_down[y][0][x];
                Ez_old_down[y][0][x] = E->z[x][y];
            }

            // calculate ABC @ y=NY
            E->z[x][NY - 1] = cf0 * (E->z[x][NY - 3] + Ez_old_up[0][1][x])
                + cf1 * (Ez_old_up[0][0][x] + Ez_old_up[2][0][x] - E->z[x][NY - 2] - Ez_old_up[1][1][x])
                + cf2 * Ez_old_up[1][0][x]
                - Ez_old_up[2][1][x];
            // update old values @ y=0
            for (y = 0; y < 3; y++) {
                Ez_old_up[y][1][x] = Ez_old_up[y][0][x];
                Ez_old_up[y][0][x] = E->z[x][NY - 1 - y];
            }
        }





//        for (int)

        // now write relevant data to file
        if (myfile.is_open())
        {
            std::cout << "Writing to file...[t=" << t << "]" << std::endl;
            for (x = 0; x < NX; x++) {
                for (y = 0; y <NY; y++) {
                    myfile << t << " " << x << " " << y << " " << E->z[x][y] << " " << params->ec[x][y] << std::endl;
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

