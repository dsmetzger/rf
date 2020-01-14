#Xk =  sum 0:N-1(xn*e^-i*2*pi*k*n/N)


import math

N=32


print "void fft_"+str(N)+"(float* x, float* X){"

layer = 0
layers = [[]]

step = 1
N1 = N
while(N1!=1):
    if (layer ==0):
        for k in range(0,N):
            layers[layer].append(k)
            #print "X["+str(k)+"]=x["+str(k)+"];"
    else:
        layers.append([])
        #even
        for k in range(0,N,2):
            layers[layer].append(layers[layer-1][k])
        #odd
        for k in range(1,N,2):
            layers[layer].append(layers[layer-1][k])
    layer+=1
    N1/=2


while(N1!=N):
    
    layer-=1
    N1*=2
print "}"
