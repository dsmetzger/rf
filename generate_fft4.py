#Xk =  sum 0:N-1(xn*e^-i*2*pi*k*n/N)


import math

N=32


print "void fft_"+str(N)+"(float* x, float* X){"

layer = 0
layers = [[]]

step = 1
while(N!=1):
    if (layer ==0):
        for k in range(0,N):
            layers[layer].append(k)
            #print "X["+str(k)+"]=x["+str(k)+"];"
    else:
        layers.append([])
        for x in range(0,step):
            for k in range(0,N,step):
                layers[layer].append(layers[layer-1][])
    layer+=1
    N/=2
    step*=2


print "}"
