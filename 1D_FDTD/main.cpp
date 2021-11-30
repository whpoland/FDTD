#include <iostream>
#include <fstream>
#include <cmath>
//#include "matplotlib-cpp-master/matplotlibcpp.h"

//namespace plt = matplotlibcpp;

/* define mathematical constants */
#define C 300000000
#define EPS_0 0.00000000000885
#define MU_0 4*M_PI*pow(10,-7)

/* CUSTOMIZABLE VALUES */
#define FREQ 300*pow(10,6)
#define EPS_R 1
#define MU_R 1

/* calculate some useful values based on inputs */
#define EPS EPS_0*EPS_R
#define MU MU_0*MU_R
#define VEL C/sqrt(EPS_R)
#define LAMBDA VEL/FREQ
#define PERIOD 1/FREQ
#define IMP sqrt(MU/EPS)
#define CRNT 1 // Courant Number

/* define bounds of our simulation */
#define DX (LAMBDA/10.)
#define DT (DX/C)
#define X_MAX 100
#define Y_MAX X_MAX
#define Z_MAX Y_MAX
#define T_MAX 200
#define N 400 // the number of elements we will track at a time
//#define N (int)ceil(X_MAX/DX)

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

// Struct to contain important info
struct data {

};

/* function headers */
double numDeriv(double* arr, int i, int delta);
double updateEx(int t, int delta, field* E, field *H, current* J, current* M);
double updateEy(int t, int delta, field* E, field *H, current* J, current* M);
double updateEz(int t, int delta, field* E, field *H, current* J, current* M);
double updateHx(int t, int delta, field* E, field *H, current* J, current* M);
double updateHy(int t, int delta, field* E, field *H, current* J, current* M);
double updateHz(int t, int delta, field* E, field *H, current* J, current* M);
void compute(field *E, field *H, current *J, current *M);

double output[N];

int main() {
//    std::cout << "Hello, World!" << std::endl;
    field E = {{0.},{0.},{0.},{0.}};
    field H = {{0.},{0.},{0.},{0.}};
    current J = {{0.},{0.},{0.}};
    current M = {{0.},{0.},{0.}};

    // run the simulation and store in file
    compute(&E, &H, &J, &M);

    // graph the data from the file
    //plot();

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

void compute(field *E, field *H, current *J, current *M) {
    // time-stepping loop
    int t;
    std::ofstream myfile ("/Users/williampoland/Library/Mobile Documents/com~apple~CloudDocs/Downloads/FDTDadd.txt");

    for (t = 0; t < T_MAX; t++) {
        // step 1 - [TE] Hx, Hy, Ez
        int x;
        for (x = 0; x < N - 1; x++)
            H->y[x] += (E->z[x+1] - E->z[x]) * DT / MU / DX;
         //never reaches Hy[N-1], so this node stays 0 -- PMC

        // step 2 - [TM] Ex, Ey, Hz
        for (x = 1; x < N; x++)
            E->z[x] += (H->y[x] - H->y[x-1]) * DT / EPS / DX;

        /* hardwire a source node */
        // this makes Ez[0] follow a  wave shape, forcing the simulation to create a wave
        //E->z[0] = exp(-(t - 30.) * (t - 30.) / 100.);

        /* additive source node */
        E->z[N - 1] += exp(-(t - 30.) * (t - 30.) / 100.);
        //    output[t] = E->z[N-1];

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
    std::cout << "Finished Writing to File" << std::endl;
}

