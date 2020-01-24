#include <math.h>
#include <complex>
#include <iostream>

#define PI       3.1415926535897932384626433832795    // PI for sine/cos calculations
#define TWOPI    6.283185307179586476925286766559     // 2*PI for sine/cos calculations
#define log10_2  0.30102999566398119521373889472449   // Log10 of 2 
#define log10_2_INV 3.3219280948873623478703194294948 // 1/Log10(2)

using namespace std;

// Returns true if N is a power of 2
bool isPwrTwo(int N, int *M)
{
    *M = (int)ceil(log10((double)N) * log10_2_INV);// M is number of stages to perform. 2^M = N
    int NN = (int)pow(2.0, *M);
    if ((NN != N) || (NN == 0)) // Check N is a power of 2. 
        return false;

    return true;
}

int get_index(const void* val, const void* base){
    //cout<<"*********"<<endl;
    //cout<<val<<endl;
    //cout<<base<<endl;
    //cout<<((unsigned long)(val)%(unsigned long)(base))/4<<endl;
    return ((unsigned long)(val)-(unsigned long)(base))/4;
    //return ((unsigned long)(val)%(unsigned long)(base))/4;
}


int get_int(float val){
    return abs(int(1.0e9*val));
}


void print_rad2FFT(int N, std::complex<float> *x, std::complex<float> *DFT)
{
    int M = 0;
    cout << "void fft_"<<N<<"(float* x, float* X){"<<endl;
    cout << "float te_re, te_im;"<<endl;
    cout << "float te_mula, te_mulb;"<<endl;
    // Check if power of two. If not, exit        
    if (!isPwrTwo(N, &M))
        throw "Rad2FFT(): N must be a power of 2 for Radix FFT";

    // Integer Variables

    int BSep;                  // BSep is memory spacing between butterflies
    int BWidth;                // BWidth is memory spacing of opposite ends of the butterfly
    int P;                     // P is number of similar Wn's to be used in that stage
    int j;                     // j is used in a loop to perform all calculations in each stage
    int stage = 1;             // stage is the stage number of the FFT. There are M stages in total (1 to M).
    int HiIndex;               // HiIndex is the index of the DFT array for the top value of each butterfly calc 
    unsigned int iaddr;        // bitmask for bit reversal 
    int ii;                    // Integer bitfield for bit reversal (Decimation in Time)
    int MM1 = M - 1;

    unsigned int i;
    int l;
    unsigned int nMax = (unsigned int)N;

    // Double Precision Variables
    double TwoPi_N = TWOPI / (double)N;    // constant to save computational time.  = 2*PI / N
    double TwoPi_NP;

    // complex Variables (See 'struct complex')
    std::complex<float> WN;               // Wn is the exponential weighting function in the form a + jb
    std::complex<float> TEMP;             // TEMP is used to save computation in the butterfly calc
    std::complex<float> *pDFT = DFT;      // Pointer to first elements in DFT array
    std::complex<float> *pLo;             // Pointer for lo / hi value of butterfly calcs
    std::complex<float> *pHi;
    std::complex<float> *pX;              // Pointer to x[n]

    void * base_x;
    void * base_X;
    base_x = x;
    base_X = DFT;

    // Decimation In Time - x[n] sample sorting
    
    for (i = 0; i < nMax; i++, DFT++)
    {
        pX = x + i;             // Calculate current x[n] from base address *x and index i.
        ii = 0;                 // Reset new address for DFT[n]
        iaddr = i;              // Copy i for manipulations
        for (l = 0; l < M; l++) // Bit reverse i and store in ii...
        {
            if (iaddr & 0x01)     // Detemine least significant bit
                ii += (1 << (MM1 - l));    // Increment ii by 2^(M-1-l) if lsb was 1
            iaddr >>= 1;                // right shift iaddr to test next bit. Use logical operations for speed increase
            if (!iaddr)
                break;
        }
        DFT = pDFT + ii;        // Calculate current DFT[n] from base address *pDFT and bit reversed index ii    
        DFT->real( pX->real());       // Update the complex array with address sorted time domain signal x[n]
        DFT->imag( pX->imag());       // NB: Imaginary is always zero
        //cout<<get_index(pX, base_x)<<endl;
        //cout<<get_index(pX+1, base_x)<<endl;
        //cout<<get_index(DFT, base_X)<<endl;
        //cout<<get_index(DFT+1, base_X)<<endl;
        cout<<"X["<<get_index(DFT, base_X)<<"]"<<"="<< "x["<< get_index(pX, base_x)<<"]" << ";"<<endl;
        cout<<"X["<<get_index(DFT, base_X)+1<<"]"<<"="<< "x["<< get_index(pX, base_x)+1<<"]" << ";"<<endl;
    }

    cout.precision(10);
    
    // FFT Computation by butterfly calculation
    for (stage = 1; stage <= M; stage++) // Loop for M stages, where 2^M = N
    {
        BSep = (int)(pow(2, stage)); // Separation between butterflies = 2^stage
        P = N / BSep;             // Similar Wn's in this stage = N/Bsep
        BWidth = BSep / 2;     // Butterfly width (spacing between opposite points) = Separation / 2.

        TwoPi_NP = TwoPi_N*P;

        for (j = 0; j < BWidth; j++) // Loop for j calculations per butterfly
        {
            if (j != 0)              // Save on calculation if R = 0, as WN^0 = (1 + j0)
            {
                //WN.real() = cos(TwoPi_NP*j)
                WN.real(  cos(TwoPi_NP*j));     // Calculate Wn (Real and Imaginary)
                WN.imag( -sin(TwoPi_NP*j));
            }
            //cout<<WN.real() <<endl;
            //cout<<WN.imag() <<endl;
            
            for (HiIndex = j; HiIndex < N; HiIndex += BSep) // Loop for HiIndex Step BSep butterflies per stage
            {
                pHi = pDFT + HiIndex;                  // Point to higher value
                pLo = pHi + BWidth;                    // Point to lower value (Note VC++ adjusts for spacing between elements)

                if (j != 0)                            // If exponential power is not zero...
                {
                    //CMult(pLo, &WN, &TEMP);          // Perform complex multiplication of Lovalue with Wn
                    TEMP.real( (pLo->real() * WN.real()) - (pLo->imag() * WN.imag()));
                    TEMP.imag( (pLo->real() * WN.imag()) + (pLo->imag() * WN.real()));
                    
                    //CSub (pHi, &TEMP, pLo);
                    pLo->real( pHi->real() - TEMP.real());       // Find new Lovalue (complex subtraction)
                    pLo->imag( pHi->imag() - TEMP.imag());
                    
                    //CAdd (pHi, &TEMP, pHi);          // Find new Hivalue (complex addition)
                    pHi->real( (pHi->real() + TEMP.real()));
                    pHi->imag( (pHi->imag() + TEMP.imag()));
                    
                    if (get_int(WN.real())==get_int(WN.imag())){
                        if (WN.imag()>0.0){
                            if (WN.real()>0.0){
                                cout<<"te_mula"<<"="<< "(X["<<get_index(pLo, base_X)<<"]"<<"*("<<abs(WN.real())<<"));" <<endl;
                                cout<<"te_mulb"<<"="<< "(X["<<get_index(pLo, base_X)+1<<"]"<<"*("<<abs(WN.real())<<"));" <<endl;
                                cout<<"te_re=te_mula-te_mulb;" <<endl;
                                cout<<"te_im=te_mula+te_mulb;" <<endl;
                            }else{
                                cout<<"te_mula"<<"="<< "(X["<<get_index(pLo, base_X)<<"]"<<"*("<<abs(WN.real())<<"));" <<endl;
                                cout<<"te_mulb"<<"="<< "(X["<<get_index(pLo, base_X)+1<<"]"<<"*("<<abs(WN.real())<<"));" <<endl;
                                cout<<"te_re=te_mulb-te_mula;" <<endl;
                                cout<<"te_im=te_mula-te_mulb;" <<endl;
                            }
                        }else{
                            if (WN.real()>0.0){
                                cout<<"te_mula"<<"="<< "(X["<<get_index(pLo, base_X)<<"]"<<"*("<<abs(WN.real())<<"));" <<endl;
                                cout<<"te_mulb"<<"="<< "(X["<<get_index(pLo, base_X)+1<<"]"<<"*("<<abs(WN.real())<<"));" <<endl;
                                cout<<"te_re=te_mula+te_mulb;" <<endl;
                                cout<<"te_im=te_mulb-te_mula;" <<endl;
                            }else{
                                cout<<"te_mula"<<"="<< "(X["<<get_index(pLo, base_X)<<"]"<<"*("<<(WN.real())<<"));" <<endl;
                                cout<<"te_mulb"<<"="<< "(X["<<get_index(pLo, base_X)+1<<"]"<<"*("<<(WN.real())<<"));" <<endl;
                                cout<<"te_re=te_mula-te_mulb;" <<endl;
                                cout<<"te_im=te_mula+te_mulb;" <<endl;
                            }
                        }
                    }else if (get_int(WN.real())==0){
                        if (WN.imag()>0.0){
                            cout<< "bad"<<endl;
                            return;
                        }
                        cout<<"te_re"<<"="<< "X["<<get_index(pLo, base_X)+1<<"];" <<endl;
                        cout<<"te_im"<<"="<< "-X["<<get_index(pLo, base_X)<<"];" <<endl;
                    }else if (get_int(WN.imag())==0){
                        if (WN.real()>0.0){
                            cout<< "bad"<<endl;
                            return;
                        }
                        cout<<"te_re"<<"="<< "-X["<<get_index(pLo, base_X)<<"];"<<endl;
                        cout<<"te_im"<<"="<< "-X["<<get_index(pLo, base_X)+1<<"];" <<endl;
                    }else{
                        cout<<"te_re"<<"="<< "(X["<<get_index(pLo, base_X)<<"]"<<"*("<<WN.real()<<"))-"<< "(X["<<get_index(pLo, base_X)+1<<"]"<<"*("<<WN.imag() << "));" <<endl;
                        cout<<"te_im"<<"="<< "(X["<<get_index(pLo, base_X)<<"]"<<"*("<<WN.imag()<<"))+"<< "(X["<<get_index(pLo, base_X)+1<<"]"<<"*("<<WN.real() << "));" <<endl;
                    }
                    
                    cout<<"X["<<get_index(pLo, base_X)<<"]"<<"="<< "X["<<get_index(pHi, base_X)<<"]"<<"-"<< "te_re" << ";" <<endl;
                    cout<<"X["<<get_index(pLo, base_X)+1<<"]"<<"="<< "X["<<get_index(pHi, base_X)+1<<"]"<<"-"<< "te_im" << ";" <<endl;
                    
                    cout<<"X["<<get_index(pHi, base_X)<<"]"<<"="<< "X["<<get_index(pHi, base_X)<<"]"<<"+"<< "te_re" << ";" <<endl;
                    cout<<"X["<<get_index(pHi, base_X)+1<<"]"<<"="<< "X["<<get_index(pHi, base_X)+1<<"]"<<"+"<< "te_im" << ";" <<endl;
                }
                else
                {
                    TEMP.real( pLo->real());
                    TEMP.imag( pLo->imag());
                    cout<<"te_re"<<"="<< "X["<<get_index(pLo, base_X)<<"];" <<endl;
                    cout<<"te_im"<<"="<< "X["<<get_index(pLo, base_X)+1<<"];" <<endl;

                    //CSub (pHi, &TEMP, pLo);
                    pLo->real( pHi->real() - TEMP.real());       // Find new Lovalue (complex subtraction)
                    pLo->imag( pHi->imag() - TEMP.imag());
                    cout<<"X["<<get_index(pLo, base_X)<<"]"<<"="<< "X["<<get_index(pHi, base_X)<<"]"<<"-"<< "te_re" << ";" <<endl;
                    cout<<"X["<<get_index(pLo, base_X)+1<<"]"<<"="<< "X["<<get_index(pHi, base_X)+1<<"]"<<"-"<< "te_im" << ";" <<endl;

                    //CAdd (pHi, &TEMP, pHi);          // Find new Hivalue (complex addition)
                    pHi->real( (pHi->real() + TEMP.real()));
                    pHi->imag( (pHi->imag() + TEMP.imag()));
                    cout<<"X["<<get_index(pHi, base_X)<<"]"<<"="<< "X["<<get_index(pHi, base_X)<<"]"<<"+"<< "te_re" << ";" <<endl;
                    cout<<"X["<<get_index(pHi, base_X)+1<<"]"<<"="<< "X["<<get_index(pHi, base_X)+1<<"]"<<"+"<< "te_im" << ";" <<endl;
                }
            }
        }
    }
    cout << "}"<<endl;
}

void rad2FFT(int N, std::complex<float> *x, std::complex<float> *DFT, int fake_sin = 0)
{
    int M = 0;

    // Check if power of two. If not, exit        
    if (!isPwrTwo(N, &M))
        throw "Rad2FFT(): N must be a power of 2 for Radix FFT";

    // Integer Variables

    int BSep;                  // BSep is memory spacing between butterflies
    int BWidth;                // BWidth is memory spacing of opposite ends of the butterfly
    int P;                     // P is number of similar Wn's to be used in that stage
    int j;                     // j is used in a loop to perform all calculations in each stage
    int stage = 1;             // stage is the stage number of the FFT. There are M stages in total (1 to M).
    int HiIndex;               // HiIndex is the index of the DFT array for the top value of each butterfly calc 
    unsigned int iaddr;        // bitmask for bit reversal 
    int ii;                    // Integer bitfield for bit reversal (Decimation in Time)
    int MM1 = M - 1;

    unsigned int i;
    int l;
    unsigned int nMax = (unsigned int)N;

    // Double Precision Variables
    double TwoPi_N = TWOPI / (double)N;    // constant to save computational time.  = 2*PI / N
    double TwoPi_NP;

    // complex Variables (See 'struct complex')
    std::complex<float> WN;               // Wn is the exponential weighting function in the form a + jb
    std::complex<float> TEMP;             // TEMP is used to save computation in the butterfly calc
    std::complex<float> *pDFT = DFT;      // Pointer to first elements in DFT array
    std::complex<float> *pLo;             // Pointer for lo / hi value of butterfly calcs
    std::complex<float> *pHi;
    std::complex<float> *pX;              // Pointer to x[n]


    // Decimation In Time - x[n] sample sorting
    for (i = 0; i < nMax; i++, DFT++)
    {
        pX = x + i;             // Calculate current x[n] from base address *x and index i.
        ii = 0;                 // Reset new address for DFT[n]
        iaddr = i;              // Copy i for manipulations
        for (l = 0; l < M; l++) // Bit reverse i and store in ii...
        {
            if (iaddr & 0x01)     // Detemine least significant bit
                ii += (1 << (MM1 - l));    // Increment ii by 2^(M-1-l) if lsb was 1
            iaddr >>= 1;                // right shift iaddr to test next bit. Use logical operations for speed increase
            if (!iaddr)
                break;
        }
        DFT = pDFT + ii;        // Calculate current DFT[n] from base address *pDFT and bit reversed index ii    
        DFT->real( pX->real());       // Update the complex array with address sorted time domain signal x[n]
        DFT->imag( pX->imag());       // NB: Imaginary is always zero
    }

    // FFT Computation by butterfly calculation
    for (stage = 1; stage <= M; stage++) // Loop for M stages, where 2^M = N
    {
        BSep = (int)(pow(2, stage)); // Separation between butterflies = 2^stage
        P = N / BSep;             // Similar Wn's in this stage = N/Bsep
        BWidth = BSep / 2;     // Butterfly width (spacing between opposite points) = Separation / 2.

        TwoPi_NP = TwoPi_N*P;

        for (j = 0; j < BWidth; j++) // Loop for j calculations per butterfly
        {   
            if (fake_sin !=1){
                if (j != 0)              // Save on calculation if R = 0, as WN^0 = (1 + j0)
                {
                    //WN.real() = cos(TwoPi_NP*j)
                    WN.real(  cos(TwoPi_NP*j));     // Calculate Wn (Real and Imaginary)
                    WN.imag( -sin(TwoPi_NP*j));
                }
                
            }else{
                WN.real( pX->real());
                WN.imag( pX->imag());
            }

            for (HiIndex = j; HiIndex < N; HiIndex += BSep) // Loop for HiIndex Step BSep butterflies per stage
            {
                pHi = pDFT + HiIndex;                  // Point to higher value
                pLo = pHi + BWidth;                    // Point to lower value (Note VC++ adjusts for spacing between elements)

                if (j != 0)                            // If exponential power is not zero...
                {
                    //CMult(pLo, &WN, &TEMP);          // Perform complex multiplication of Lovalue with Wn
                    TEMP.real( (pLo->real() * WN.real()) - (pLo->imag() * WN.imag()));
                    TEMP.imag( (pLo->real() * WN.imag()) + (pLo->imag() * WN.real()));

                    //CSub (pHi, &TEMP, pLo);
                    pLo->real( pHi->real() - TEMP.real());       // Find new Lovalue (complex subtraction)
                    pLo->imag( pHi->imag() - TEMP.imag());

                    //CAdd (pHi, &TEMP, pHi);          // Find new Hivalue (complex addition)
                    pHi->real( (pHi->real() + TEMP.real()));
                    pHi->imag( (pHi->imag() + TEMP.imag()));
                }
                else
                {
                    TEMP.real( pLo->real());
                    TEMP.imag( pLo->imag());

                    //CSub (pHi, &TEMP, pLo);
                    pLo->real( pHi->real() - TEMP.real());       // Find new Lovalue (complex subtraction)
                    pLo->imag( pHi->imag() - TEMP.imag());

                    //CAdd (pHi, &TEMP, pHi);          // Find new Hivalue (complex addition)
                    pHi->real( (pHi->real() + TEMP.real()));
                    pHi->imag( (pHi->imag() + TEMP.imag()));
                }
            }
        }
    }
}
