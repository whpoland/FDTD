#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
//#include "matplotlib-cpp-master/matplotlibcpp.h"

//namespace plt = matplotlibcpp;

/* define mathematical constants */
#define C 3e8
#define EPS_0 8.85e-12
#define MU_0 (M_PI*4e-7)

/* CUSTOMIZABLE VALUES */
#define FREQ 5.9e9 // used for scaling simulation parameters (DT, DX, etc.)
#define PROP_FREQ FREQ // used for source wave

#define FRAC 10. // fraction of wavelength for spatial step (spatial resolution of grid) *Should be 20 for high accuracy*
#define CNO (0.7071) // Courant Number - stability factor


/* calculate useful values */
#define LAMBDA (C/FREQ) // wavelength
#define OMEGA (2*M_PI*FREQ)
#define PROP_OMEGA (2*M_PI*PROP_FREQ)
//#define BETA (OMEGA/C) // phase constant

/* define bounds of our simulation */
#define DX (LAMBDA/FRAC)
#define DY DX
#define DT (DX/C*CNO)
//#define X_MAX 100
//#define Y_MAX X_MAX
//#define Z_MAX Y_MAX
//#define N (int)ceil(X_MAX/DX)


#define T_MAX 1000 // how long the simulation will run
//#define T_MAX 20 * round(1.0/FREQ/DT) //in terms of number of periods
//#define N (int)(FRAC * 3 / 2) // spatial size of grid
#define N 1000
#define NX N
#define NY NX

#define GEO_PRECISION double
#define PRECISION float

//const int T_MAX, N, NX, NY;

// Create a struct template for E,H fields (to contain x,y,z, components)
struct field {
    std::array<std::array<PRECISION,NX>,NY> x;
    std::array<std::array<PRECISION,NX>,NY> y;
    std::array<std::array<PRECISION,NX>,NY> z;
};

// Parameters of the simulation
struct parameters {
    std::array<std::array<GEO_PRECISION,NX>,NY> eps; // relative permittivity of medium
    std::array<std::array<GEO_PRECISION,NX>,NY> mu; // relative permeability of medium
    std::array<std::array<GEO_PRECISION,NX>,NY> ec; // electrical conductivity of medium
    std::array<std::array<GEO_PRECISION,NX>,NY> mc;  // magnetic conductivity of medium

    std::string test_name; // name of the output .txt file
};

/* function headers */
void compute(field *E, field *H, field *J, field *M, parameters *params);

int main() {
    std::cout << "Initializing Simulation..." << std::endl;
    std::cout << "NX: " << NX << " | NY: " << NY << " | TMAX: " << T_MAX;
    std::cout << " | DX: " << DX << " | DY: " << DY << " | DT: " << DT << std::endl;

//  Initialize fields
    field *E = new field;
    field *H = new field;
    field *J = nullptr;
    field *M = nullptr;
    parameters *params = new parameters;
    // adjust name of simulation output file
    params->test_name = "indoor2";

    int i, j; // i is x location, j is y location
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            // defaults
            E->x[i][j] = 0;
            E->y[i][j] = 0;
            E->z[i][j] = 0;
            H->x[i][j] = 0;
            H->y[i][j] = 0;
            H->z[i][j] = 0;
            params->eps[i][j] = 1;
            params->mu[i][j] = 1;
            params->ec[i][j] = 0;
            params->mc[i][j] = 0;
            // add in conditions to modify material properties at specific points (i.e. add dielectrics and conductors)


            // ECE 4124 Final: Indoor Propagation
            double wall_eps = 2.8;
            double wall_cond = 0.0115;
            double door_eps = 2.13;
            double door_cond = 0.7;
            int thickness = 3;
            int doorlen = 12;

            // indoor rooms
            if (((i >= 0 && i <= 0 + thickness) || (i >= NX-1-thickness && i <= NX-1)) || // vertical outer walls
                ((j >= 0 && j <= 0 + thickness) || (j >= NY-1-thickness && j <= NY-1)) || // horizontal outer walls
                ((j >= 2*NY/5-thickness && j <= 2*NY/5) || (j >= 3*NY/5-thickness && j <= 3*NY/5)) || // horizontal inner walls
                ((j <= 2*NY/5-thickness || j >= 3*NY/5) && ((i >= NX/3-thickness && i <= NX/3) ||   // -- vertical inner walls
                                                            (i >= 2*NX/3-thickness && i <= 2*NX/3)))  // -- vertical inner walls
                ) {
                params->eps[i][j] = wall_eps;
                params->ec[i][j] = wall_cond;
            }
            // doors
            if (((j >= 2*NY/5-thickness && j <= 2*NY/5) || (j >= 3*NY/5-thickness && j <= 3*NY/5)) &&
                    ((i >= NX/6-doorlen/2  && i <= NX/6+doorlen/2) ||
                     (i >= 3*NX/6-doorlen/2 && i <= 3*NX/6+doorlen/2) ||
                     (i >= 5*NX/6-doorlen/2 && i <= 5*NX/6+doorlen/2))
                     ) {
                params->eps[i][j] = door_eps;
                params->ec[i][j] = door_cond;
            }


            //***************************** properties for ECE 4124 Project 2 ******************************
            // diffraction pattern
//            if (i == NX / 2) {
//                if ((j < NX / 3 + 5 && j > NX / 3 - 5) || (j < 2 * NX / 3 + 5 && j > 2 * NX / 3 - 5)) {
//                    continue;
//                }
//                else {
//                    params.ec[i][j] = 1e9;
//                }
//            }
//            if (i == NX / 2 && (j/20) % 2 == 0) {
//                    params.ec[i][j] = 1e9;
//            }

            // knife edge
//            if (i == NX / 2 && j < NY / 2) { // half plane
////            if (i == NX / 2 && j < 3 * NY / 4) { // 3/4 plane
//                params.ec[i][j] = 1e9;
//            }
            // round edge
//            if ((i-NX/2)*(i-NX/2) + (j-NY/2)*(j-NY/2) <= 25) {
//                params.ec[i][j] = 1e9;
//            }
//            if (j < NY/2 && (i == NX/2-5 || i == NX/2+5)) {
//                params.ec[i][j] = 1e9;
//            }


            // ground plane -- might have issues
//            if (j == 1) {
//                params.ec[i][j] = 1e9;
//            }
            //*******************************************************************************************************

            // ****************** setting conductor boundaries for waveguide simulation *********************************************
//            // rectangular waveguide - define conductive boundaries
//            int x1 = params.src_x - + round(LAMBDA / 4 / DY); // starting x location
////            x1 = params.src_x - 2;
////            int x2 = NX - 1; // ending x location
////            int wav_wid = 1.4 * LAMBDA / 2 / DY; // make width >lambda/2 to avoid cutoff
//            int wav_wid = 22.86e-3 / DY; // width of WR90 waveguide at 10GHz frequency
//            int y1 = params.src_y - floor(wav_wid/2); // min y value (bottom)
//            int y2 = params.src_y + ceil(wav_wid/2); // max y value (top)
//            if ((i > floor(x1 - 1) && i < floor(x1 + 1) && j > floor(y1 - 1) && j < floor(y2 + 1))
////                    || (i > floor(x2 - 1) && i < floor(x2 + 1) && j > floor(y1 - 1) && j < floor(y2 + 1))
//                    || (j > floor(y1 - 1) && j < floor(y1 + 1) && i > floor(x1 - 1))
//                    || (j > floor(y2 - 1) && j < floor(y2 + 1) && i > floor(x1 - 1))
//                    ){
//                params.ec[i][j] = 1e9;
//            }
//            // rectangular waveguide - set dielectric thickness within metal
//            int x3 = round(NX/2);
//            int sample_thickness = 1.1*round(LAMBDA/2/DY); // multiply by factor of 1.1
//
//            if (i > x3 && i < x3 + sample_thickness && j > y1 - 1 && j < y2 + 1)
////                params.eps[i][j] = 10;
//                ;
            // ************************************************************************************************************
        }
    }

    // run the simulation and store in file
    compute(E, H, J, M, params); // J and M are currently declared as nullptr

    return 0;
}

void compute(field *E, field *H, field *J, field *M, parameters *params) {
    // time-stepping loop
    std::cout << "Starting Simulation..." << std::endl;
    int t;
    std::string file_path = "/Users/williampoland/FDTD/2D/raw/";
    std::string output_file_name = "/Users/williampoland/FDTD/2D/raw/" + params->test_name + ".txt";
    std::string geometry_file_name = file_path + params->test_name + "_geometry.txt";
    std::ofstream myfile (output_file_name); // file to output E field data
    std::ofstream geofile (geometry_file_name); // file to output geometry and simulation parameters

    geofile << "FREQ T_MAX NX NY DT DX DY" << std::endl;
    geofile << PROP_FREQ << " " << T_MAX << " " << NX << " " << NY << " " << DT << " " << DX << " " << DY << std::endl;
    geofile << "(eps_r, ele_cond)" << std::endl;

    myfile << "(Ex, Ey, Ez)" << std::endl;

    // FDTD simulation
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
                C_Ey_H = DT / EPS_0 / params->eps[x][y] / DY / (1 + cf_Ey);

                E->y[x][y] = C_Ey_E * E->y[x][y] + C_Ey_H * (H->z[x][y] - H->z[x-1][y]);
            }
        }

        // Update Hz
        double cf_Hz, C_Hz_E, C_Hz_H;
        for (x = 0; x < NX - 1; x++) {
            for (y = 0; y < NY - 1; y++) {
                cf_Hz = params->mc[x][y] * DT / 2 / MU_0 / params->mu[x][y];
                C_Hz_H = (1 - cf_Hz) / (1 + cf_Hz);
                C_Hz_E = DT / MU_0 / params->mu[x][y] / DX / (1 + cf_Hz);

                H->z[x][y] = C_Hz_H * H->z[x][y] + C_Hz_E * ((E->x[x][y+1] - E->x[x][y]) - (E->y[x+1][y] - E->y[x][y]));
            }
        }
        // ********************* SET SOURCES ************************************************************************************
        /* create source node */
        // Note: source nodes cannot be on boundaries x=0,x=NX,y=0,y=NY because of the absorbing boundaries
        // Source nodes for grid edges should be placed one node from the edge of grid (i.e. x=1 or y=NY-1)
        // Can create hardwired or additive sources: hardwired - enforce field at a point; additive - add to field at a point
        int delay = 0;
        int width = N/16;
        int period = 1 / PROP_FREQ / DT;
        // hard-wired source

        // indoor propagation simulation (TMz & TEz)
        float src_mag = 0.632;
        E->z[NX/2][NY/2] += src_mag * sin(PROP_OMEGA * t * DT);
        E->x[NX/2][NY/2] += E->z[NX/2][NY/2];
        E->y[NX/2][NY/2] += E->z[NX/2][NY/2];

//        E->z[1][1] = exp(-(t - 30.) * (t - 30.) / 100.);
//        E->z[NX/2][NY/2] = 1.5 * exp(-(t - delay) * (t - delay) / width);
        // source for waveguide sim
//        E->z[params->src_x][params->src_y] = sin(PROP_OMEGA*t*DT);
        // vertical plane wave
//        for (int i = 0; i < NX; i++) {
////            E->z[i][0] = sin(t);
//            E->z[i][1] = exp(-(t - delay) * (t - delay) / width);
//        }
        // spherical wave
//        E->z[1][3*NY/4] = sin(PROP_OMEGA * t * DT);
//         *** horizontal plane wave ***
//        for (int j = 1; j < NY; j++) {
//            E->z[1][j] = sin(PROP_OMEGA * t * DT);
            // turn off source after one period
//            if (t <= period) {
//                E->z[1][j] = sin(PROP_OMEGA * t * DT);
//            }
//            else {
//                E->z[1][j] = 0;
//            }
//           E->z[1][j] = exp(-(t - delay) * (t - delay) / width);
//        }
        // additive source
//        E->z[NX/2][NY/2] += 1.5*exp(-(t - delay) * (t - delay) / width);

        // ***************************************************************************************************

        // ****************************************   ABCs    ************************************************
        // initialize variables for boundary conditions
        // these will store the past values of for edges of grid
        // 1st index - displacement from edge of boundary
        // 2nd index - previous time displacement
        // 3rd index - location along boundary
        const int num_disp = 3, time_disp = 2;
        double
//        Ex_old_up[num_disp][time_disp][NX] = {0.,0.}, Ex_old_down[num_disp][time_disp][NX] = {0.,0.},
//                Ex_old_left[num_disp][time_disp][NY] = {0.,0.}, Ex_old_right[num_disp][time_disp][NY] ={0.,0.},
//                Ey_old_up[num_disp][time_disp][NX] = {0.,0.}, Ey_old_down[num_disp][time_disp][NX] = {0.,0.},
//                Ey_old_left[num_disp][time_disp][NY] = {0.,0.}, Ey_old_right[num_disp][time_disp][NY] ={0.,0.},
                Ez_old_up[num_disp][time_disp][NX] = {0.,0.}, Ez_old_down[num_disp][time_disp][NX] = {0.,0.},
                Ez_old_left[num_disp][time_disp][NY] = {0.,0.}, Ez_old_right[num_disp][time_disp][NY] ={0.,0.};

        /* second order ABCs for Ex,Ey,Ez @ x = 0, x = NX, y = 0, y = NY */

        // calculate coefficients
        cf_Ez = params->ec[0][0] * DT / 2 / EPS_0 / params->eps[0][0];
        C_Ez_H = DT / EPS_0 / params->eps[0][0] / DY / (1 + cf_Ez);

        cf_Hy = params->mc[0][0] * DT / 2 / MU_0 / params->mu[0][0];
        C_Hy_E = DT / MU_0 / params->mu[0][0] / DY / (1 + cf_Hy);

        double temp1_TM = sqrt(C_Ez_H * C_Hy_E);
        double temp2_TM = 1. / temp1_TM + 2. + temp1_TM;
        double cf0 = -(1. / temp1_TM - 2. + temp1_TM) / temp2_TM;
        double cf1 = -2. * (temp1_TM - 1. / temp1_TM) / temp2_TM;
        double cf2 = 4. * (temp1_TM + 1. / temp1_TM) / temp2_TM;


        // ABCs @ x=0 & x=NX
        for (y = 0; y < NY; y++) {
            // calculate ABC @ x=0
//            E->y[0][y] = 0;

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
//            E->x[NX - 1][y] = 0;

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
//            E->x[x][0] = 0;

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
//            E->x[x][NY - 1] = 0;

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






        // now write relevant data to file
        if (t == 0 && geofile.is_open()) { // only write our geometry details once at time t=0
            std::cout << "Writing to geometry file" << std::endl;
            for (x = 0; x < NX; x++) {
                for (y = 0; y <NY; y++) {
                    geofile << params->eps[x][y] << " " << params->ec[x][y] << std::endl;
                }
            }
            geofile.close(); // close file when done
            std::cout << "Finished writing to geometry file" << std::endl;
        }
        // write field data to file for each time t
        if (myfile.is_open()) {
            std::cout << "Writing to file...[t=" << t << "]" << std::endl;
            for (x = 0; x < NX; x++) {
                for (y = 0; y <NY; y++) {
//                    myfile << t << " " << x << " " << y << " " << E->z[x][y] << " " << params->eps[x][y] << " " << params->ec[x][y] << std::endl;
                    myfile << E->x[x][y] << " " << E->y[x][y] << " " << E->z[x][y]  << std::endl;

//                    myfile << t << " " << x << " " << y << " "
//                            << E->z[x][y] << " "  << sqrt((E->x[x][y]*E->x[x][y]) + (E->y[x][y]*E->y[x][y])) << " "
//                            << params->eps[x][y] << " " << params->ec[x][y] << std::endl;

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
    delete E;
    delete H;
    delete params;
    std::cout << "Finished Writing to File: " << output_file_name << std::endl;
    std::cout << "Total entries: " << NX * NY * T_MAX << std::endl;
    std::cout << "Simulation Complete." << std::endl;
}

