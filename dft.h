/*
Library for discrete fourier transform algorithims. Everything should be either templated or have functions for complex floats/shorts.

*/

#include <vector>


float cas(float m_rad){
    const float c1 = 1.41421356237;
    const float c2 = 0.78539816339;
    return c1*cos(m_rad-c2);
}

/*
//class LUT
class casLUT{
private:

public:

}
*/

class DHT{
private:
    std::vector<float> casLUT;
    int N
public:
    DHT(int N):N(N){
        casLUT.resize(N*N);
        int idx = 0;
        double factor = 2.0*M_PI/double(N);
        double factor1;
        for (int k=0; k<N; ++k){
            factor1 = factor*k;
            for (int n=0; n<N; ++n){
                casLUT[idx]= cas(factor1*n);
                ++idx;
            }
        }
    }
    
    void brute_dht(float* x, float* Xk){
        int idx = 0;
        for (int k=0; k<N; ++k){
            Xk[k] = 0;
            for (int n=0; n<N; ++n){
                Xk[k]+= x[n]*casLUT[idx];
                ++idx;
            }
        }
    }
    
    void convolve_dht(){
        
    }
}




