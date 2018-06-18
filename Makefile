INCLUDE_PATH= ./eigen-3.3.4
CFLAGS=$(foreach d, $(INCLUDE_PATH), -I$d) -Wall -std=c++11 -O3
VRPN_LIBS = -lvrpn -lquat -pthread
OPTI_FLAG = -O3
GTK_LIBS = -rdynamic `pkg-config --cflags gtk+-3.0`  `pkg-config --libs gtk+-3.0` 

test_FFTSO3: test_FFTSO3.o fdcl_FFTSO3.o fdcl_FFTSO3_matrix.o misc_matrix_func.o 
	g++ -o test_FFTSO3  test_FFTSO3.o fdcl_FFTSO3.o fdcl_FFTSO3_matrix.o misc_matrix_func.o -lm $(CFLAGS) $(OPTI_FLAGS)

test_FFTSO3.o: test_FFTSO3.cpp
	g++ -c test_FFTSO3.cpp $(CFLAGS) $(OPTI_FLAGS)

fdcl_FFTSO3.o: fdcl_FFTSO3.cpp fdcl_FFTSO3.hpp
	g++ -c fdcl_FFTSO3.cpp $(CFLAGS) $(OPTI_FLAGS)
	
fdcl_FFTSO3_matrix.o: fdcl_FFTSO3_matrix.hpp fdcl_FFTSO3_matrix.cpp
	g++ -c fdcl_FFTSO3_matrix.cpp $(CFLAGS) $(OPTI_FLAGS)

misc_matrix_func.o: misc_matrix_func.cpp misc_matrix_func.h
	g++ -c misc_matrix_func.cpp $(CFLAGS) $(OPTI_FLAGS)

