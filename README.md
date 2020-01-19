# rf



g++ test_dft.cc -Wall -O0 -lm -std=c++11
10*4096/903931 nanoseconds
g++ test_dft.cc -Wall -O3 -lm -std=c++11 -pthread -lpthread -pipe
10*4096/462123 nanoseconds


g++ test_dft_32768.cc -o 32768.out -Wall -O0 -lm -std=c++11 -pthread -lpthread -pipe
g++ test_dft_32768.cc -o 32768.out -Wall -O3 -lm -std=c++11 -pthread -lpthread -pipe
