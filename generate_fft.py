#Xk =  sum 0:N-1(xn*e^-i*2*pi*k*n/N)

'''
example

in_real = xn[n];++n;in_imag = xn[n];++n;

out[20] = in_real;
re = 1.1;
im = .5;
out[20] += re*in_real - im*in_imag;
re = .1;
im = .2;
out[20] += re*in_real - im*in_imag;
++k;
out[k] = re*in_imag + im*in_real;
++k;

'''


import math

N=4096


#section off by k*n%N and (xn*e^-i*2*pi*k*n/N)
exp_strs = {}

for k in range(0,N):
    for n in range(0,N):
        exp_strs[str(n)+"_"+str(k*n%N)]=""

factor = 2.0*math.pi/float(N);
factor1 = 0.0;
for k in range 0:N:
    factor1 = float(k)*factor
    for n in range 0:N:
        out_ix = (n*k)%N
        key_to_str = str(n)+"_"+str(k*n%N)
        exp_strs[key_to_str] += "out["+str(2*out_ix)+"] = out_real;out["+str(2*out_ix+1)+"] = out_imag;"


for k in range 0:N:
    print "in_real = xn[n];++n;in_imag = xn[n];++n;"
    for n in range(0,N):
        key_to_str = str(n)+"_"+str(k*n%N)
        print exp_strs[key_to_str]
    




print "void fft_"+str(N)+"(){"
print "int n = 0;int k = 0;float re,im;float in_real,in_imag;float out_real,out_imag;"

print "}"
