/*
Library for discrete fourier transform algorithims. Everything should be either templated or have functions for complex floats/shorts.

*/

#include <vector>


void cas(float m_rad){
    
}

//class LUT
class casLUT{
private:

public:

}


//DHT functions
void brute_dht(short* x, float* Xk, const int N){
    double factor = 2.0*M_PI/double(N);
    double factor1;
    for (int k=0; k<N; ++k){
        factor1 = factor*k;
        Xk[k] = 0;
        for (int n=0; n<N; ++n){
            Xk[k]+= x[n]*cas(factor1*n)
        }
    }
}

void convolve_dht(){

}
