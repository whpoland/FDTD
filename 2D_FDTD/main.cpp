#include <iostream>
#include <fstream>
#include <cmath>

/* define mathematical constants */
#define C 3*pow(10,8)
#define EPS_0 8.85*pow(10,-12)
#define MU_0 4*M_PI*pow(10,-7)

/* CUSTOMIZABLE VALUES */
#define N 10 // the number of elements we will track at a time
#define FREQ 300*10^6
#define EPS_R 1
#define MU_R 1
#define NUM_LENGTHS 50 // how many wavelengths
#define ELEC_COND 0
#define MAG_COND 0

/* calculate some useful values based on inputs */
#define EPS EPS_0*EPS_R
#define MU MU_0*MU_R
#define VEL C/sqrt(E_R)
#define LAMBDA VEL/FREQ
#define PERIOD 1/FREQ

/* define bounds of our simulation */
#define X_MAX NUM_LENGTHS*LAMBDA
#define Y_MAX X_MAX
#define Z_MAX Y_MAX
#define T_MAX NUM_LENGTHS*PERIOD*1.5

// Create a struct template for E,H,J,M fields (to contain x,y,z components)
struct field {
    double x[N];
    double y[N];
    double z[N];
};

/* function headers */
double numDeriv(double* arr, int i, int delta);
double updateEx(int t, int delta, field* E, field *H, field* J, field* M);
double updateEy(int t, int delta, field* E, field *H, field* J, field* M);
double updateEz(int t, int delta, field* E, field *H, field* J, field* M);
double updateHx(int t, int delta, field* E, field *H, field* J, field* M);
double updateHy(int t, int delta, field* E, field *H, field* J, field* M);
double updateHz(int t, int delta, field* E, field *H, field* J, field* M);
void compute(field *E, field *H, field *J, field *M);
void plot();

int main() {
    std::cout << "Hello, World!" << std::endl;
    field E;
    field H;
    field J;
    field M;

    // run the simulation and store in file
    compute(&E, &H, &J, &M);

    // graph the data from the file
    plot();

    return 0;
}

double numDeriv(double* arr, int i, int delta) {
    return (arr[i + delta/2] + arr[i - delta/2])/delta;
}

double updateEx(int t, int delta, field* E, field *H, field* J, field* M) {

}

double updateHx(int t, int delta, field* E, field *H, field* J, field* M) {
    H->x[t - delta/2] + delta/MU*(numDeriv(E->z, t, delta) + M->x[t] + MAG_COND*H->x[t])
}

void compute(field *E, field *H, field *J, field *M) {
// loop
    // step 1 - [TE] Hx, Hy, Ez
    // step 2 - [TM] Ex, Ey, Hz
};

void plot() {
    ;
}