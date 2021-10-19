#include <iostream>

#define N 10 // the number of elements we will track at a time

// Create a struct template for E and H fields (to contain x,y,z components)
struct field {
    double x[N];
    double y[N];
    double z[N];
};
int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
