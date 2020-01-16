#include <iostream>
#include <vector>
#include <complex>
#include "dft_32.h"
#include <math.h>

#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

void brute_dft(std::complex<float>* x, std::complex<float>* Xk, int N){
	int idx = 0;
	float factor = 2.0*M_PI/float(N);
	float factor1; 
	for (int k=0; k<N; ++k){
	    Xk[k] = 0;
		factor1 = factor*k;
	    for (int n=0; n<N; ++n){
			Xk[k] += x[n]*std::complex<float>(cos(factor1*n), -sin(factor1*n));
	    }
	}
}

inline std::complex<float> thexp(float factor){
	return (std::complex<float>(cos(factor), -sin(factor)));
}

void ditfft2(std::complex<float>* x, std::complex<float>* X, int N, int s){
    if (N == 2){
std::cout<< "Xe "<<&X[0] << std::endl;
std::cout<< "Xo "<<&X[s] << std::endl;
        X[0] = x[0]+x[s];
        X[s] = x[0]-x[s];
    }else{
std::cout<< "N1 "<<N << std::endl;
        std::complex<float> t;
        ditfft2(x, X, N/2, 2*s);
        ditfft2(x+s, X+s, N/2, 2*s);
std::cout<< "N2 "<<N << std::endl;
std::cout<< "s "<<s << std::endl;
        for (int k = 0; k<N/2; k+=1){
//std::cout<< "k "<<k << std::endl;
std::cout<< "Xe "<<&X[2*k*s] << std::endl;
std::cout<< "Xo "<<&X[2*k*s+s] << std::endl;
            t = X[2*k*s];
            X[k*s] = t + thexp(2*M_PI* k/N)*X[2*k*s+s];
            X[2*k*s+s*N/2] = t - thexp(2*M_PI* k/N)*X[2*k*s+s];
        }
    }
}

int main(){
	int N = 4;
	std::vector<std::complex<float> > input_vec(N,1.0);
	input_vec[N/2] = float(N);
	std::vector<std::complex<float> > o0(N,1.0);
	std::vector<std::complex<float> > o1(N,1.0);

	auto start_time = Clock::now();
	for (int x=0; x<1; ++x){
		brute_dft(&input_vec[0], &o0[0], N);
	}
	auto end_time = Clock::now();
	std::cout << "Time difference:"
      << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;

	for (int x=0; x<N; ++x){
		std::cout<<o0[x]<<",";
	}
	std::cout<<std::endl;std::cout<<std::endl;

std::cout<< "X0 "<<&o1[0] << std::endl;
std::cout<< "X1 "<<&o1[1] << std::endl;


	start_time = Clock::now();
	for (int x=0; x<1; ++x){
		ditfft2(&input_vec[0], &o1[0], N, 1);
	}
	end_time = Clock::now();
	std::cout << "Time difference:"
      << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
/*
    start_time = Clock::now();
	for (int x=0; x<10; ++x){
		fft_32((float *)&input_vec[0], (float *)&o1[0]);
	}
	end_time = Clock::now();
	std::cout << "Time difference:"
      << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
*/
	for (int x=0; x<N; ++x){
		std::cout<<o1[x]<<",";
	}
	std::cout<<std::endl;

	return 0;
}
