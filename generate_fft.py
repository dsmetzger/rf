#Xk =  sum 0:N-1(xn*e^-i*2*pi*k*n/N)

'''

example

ir = xn[n];++n;ii = xn[n];++n;

o[20] = ir;
re = 1.1;
im = .5;
o[20] += re*ir - im*ii;
re = .1;
im = .2;
o[20] += re*ir - im*ii;
++k;
o[k] = re*ii + im*ir;
++k;

'''


import math

N=32

def round(infloat):
    return int(infloat*10.0e6)

#section off by k*n%N and (xn*e^-i*2*pi*k*n/N)
exp_strs = {}
re_strs = {}
im_strs = {}
dict_phi = {}

#print "setting up dictionary"
factor = 2.0*math.pi/float(N);
factor1 = 0.0;
for k in range(0,N):
    factor1 = float(k)*factor
    for n in range(0,N):
        exp_strs[str(n)+"_"+str((k*n)%N)]="\t\t"
        dict_phi[str(n)+"_"]=[]
        
        #key_to_str_re = str(n)+"_"+str(round(math.cos(n*factor1)))
        #re_strs[key_to_str_re]="\t\t"
        #key_to_str_im = str(n)+"_"+str(round(-math.sin(n*factor1)))
        #im_strs[key_to_str_im]="\t\t"


#print "adding strings to dictionary"
factor = 2.0*math.pi/float(N);
factor1 = 0.0;
for k in range(0,N):
    #print "k ",k
    factor1 = float(k)*factor
    for n in range(0,N):
        out_ix = k
        key_to_str = str(n)+"_"+str((k*n)%N)
        #key_to_str_re = str(n)+"_"+str(round(math.cos(n*factor1)))
        #key_to_str_im = str(n)+"_"+str(round(-math.sin(n*factor1)))
        #print key_to_str
        if n==0:
            exp_strs[key_to_str] += "o["+str(2*out_ix)+"]=ore;o["+str(2*out_ix+1)+"]=oi;"
            #re_strs[key_to_str_re] += "Xk["+str(2*out_ix)+"]=ore;"
            #im_strs[key_to_str_im] += "Xk["+str(2*out_ix+1)+"]=oi;"
        else:
            exp_strs[key_to_str] += "o["+str(2*out_ix)+"]+=ore;o["+str(2*out_ix+1)+"]+=oi;"
            #re_strs[key_to_str_re] += "Xk["+str(2*out_ix)+"]+=ore;"
            #im_strs[key_to_str_im] += "Xk["+str(2*out_ix+1)+"]+=oi;"


print "#include <iostream>"
print "#include <math.h>"
print "\n"
print "\n"
print "void fft_"+str(N)+"(float* xn, float* o, float* buffer = 0){"
print "\tint n = 0;int k = 0;float re,im;float ir,ii;float ore,oi;"

print "\n"
print "\n"

#print exp_strs["0_0"]
#print "\n"
#print exp_strs["1_0"]

#quit()

print "\n"

factor = 2.0*math.pi/float(N);
factor1 = 0.0;
a_str = ""
b_str = ""
for n in range(0,N):
    factor1 = float(n)*factor
    print "\t//n, ",n
    print "\tir = xn[n];++n;ii = xn[n];++n;"
    tmp_dict = {}
    for k in range(0,N):
        #print "\t\t//k, ",k
        #print "\t\t//n*k%N, ",(k*n)%N
        '''
        key_to_str = str(n)+"_"+str((k*n)%N)
        if (exp_strs[key_to_str] != ""):
            real = math.cos(k*factor1)
            imag = -math.sin(k*factor1)
            key_to_str_re = round(real)
            key_to_str_im = round(imag)
            if key_to_str_re in tmp_dict.keys():
                a_str = tmp_dict[key_to_str_re]

            else:
                #a_str = "buffer["+str(2*k)+"]"
                #print a_str+"="+str(real)+"*ir"
                print "\t\tore="+str(real)+"*ir-"+str(imag)+"*ii;"
                print "\t\toi="+str(real)+"*ii+"+str(imag)+"*ir;"

                #print "ore = re*ir - im*ii"
                #print "oi = re*ii + im*ir"
                print exp_strs[key_to_str]
                exp_strs[key_to_str] = ""

            if key_to_str_im in tmp_dict.keys():
                
            else:
                
                print "\t\tore="+str(real)+"*ir-"+str(imag)+"*ii;"
                print "\t\toi="+str(real)+"*ii+"+str(imag)+"*ir;"

                #print "ore = re*ir - im*ii"
                #print "oi = re*ii + im*ir"
                print exp_strs[key_to_str]
                exp_strs[key_to_str] = ""
            
            print "\t\tore="+a_str+"-"+b_str+";"
            print "\t\toi="+str(real)+"*ii+"+str(imag)+"*ir;"

        '''        
        key_to_str = str(n)+"_"+str((k*n)%N)
        if (exp_strs[key_to_str] != ""):
            real = math.cos(k*factor1)
            imag = -math.sin(k*factor1)
            print "\t\tore = "+str(real)+"*ir - "+str(imag)+"*ii;"
            print "\t\toi = "+str(real)+"*ii + "+str(imag)+"*ir;"

            #print "ore = re*ir - im*ii"
            #print "oi = re*ii + im*ir"
            print exp_strs[key_to_str]
            exp_strs[key_to_str] = ""
        
        '''
        real = math.cos(k*factor1)
        imag = -math.sin(k*factor1)
        key_to_str_re = str(n)+"_"+str(round(real))
        key_to_str_im = str(n)+"_"+str(round(imag))
        if (re_strs[key_to_str_re] != ""):
            print "\t\tore = "+str(real)+"*ir - "+str(imag)+"*ii;"

            print re_strs[key_to_str_re]
            re_strs[key_to_str_re] = ""

        if (im_strs[key_to_str_im] != ""):
            print "\t\toi = "+str(real)+"*ii + "+str(imag)+"*ir;"

            print im_strs[key_to_str_im]
            im_strs[key_to_str_im] = ""
        '''
print "}"
