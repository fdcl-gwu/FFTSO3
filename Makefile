INCLUDE_PATH= /usr/local/include/eigen3
CFLAGS=$(foreach d, $(INCLUDE_PATH), -I$d) -Wall -std=c++11 -O3
VRPN_LIBS = -lvrpn -lquat -pthread
OPTI_FLAG = -O3
GTK_LIBS = -rdynamic `pkg-config --cflags gtk+-3.0`  `pkg-config --libs gtk+-3.0` 

test_FFTSO3: test_FFTSO3.o fdcl_FFTS2.o fdcl_Clebsch_Gordon_matrix.o fdcl_FFTSO3.o  fdcl_FFTSO3_matrix.o fdcl_tictoc.o misc_matrix_func.o fdcl_FFTS2_matrix.o
	g++ -o test_FFTSO3  test_FFTSO3.o fdcl_FFTS2.o fdcl_Clebsch_Gordon_matrix.o fdcl_FFTSO3.o fdcl_FFTSO3_matrix.o fdcl_tictoc.o misc_matrix_func.o fdcl_FFTS2_matrix.o -lm $(CFLAGS) $(OPTI_FLAGS)

test_FFTSO3.o: test_FFTSO3.cpp
	g++ -c test_FFTSO3.cpp $(CFLAGS) $(OPTI_FLAGS)

fdcl_FFTSO3.o: fdcl_FFTSO3.cpp fdcl_FFTSO3.hpp
	g++ -c fdcl_FFTSO3.cpp $(CFLAGS) $(OPTI_FLAGS)

fdcl_FFTS2.o: fdcl_FFTS2.cpp fdcl_FFTS2.hpp
	g++ -c fdcl_FFTS2.cpp $(CFLAGS) $(OPTI_FLAGS)

fdcl_Clebsch_Gordon_matrix.o: fdcl_Clebsch_Gordon_matrix.hpp fdcl_Clebsch_Gordon_matrix.cpp
	g++ -c fdcl_Clebsch_Gordon_matrix.cpp $(CFLAGS) $(OPTI_FLAGS)

fdcl_tictoc.o: fdcl_tictoc.cpp fdcl_tictoc.hpp
	g++ -c fdcl_tictoc.cpp $(CFLAGS) $(OPTI_FLAGS)
	
fdcl_FFTSO3_matrix.o: fdcl_FFTSO3_matrix.hpp fdcl_FFTSO3_matrix.cpp
	g++ -c fdcl_FFTSO3_matrix.cpp $(CFLAGS) $(OPTI_FLAGS)

fdcl_FFTS2_matrix.o: fdcl_FFTS2_matrix.hpp fdcl_FFTS2_matrix.cpp
	g++ -c fdcl_FFTS2_matrix.cpp $(CFLAGS) $(OPTI_FLAGS)

misc_matrix_func.o: misc_matrix_func.cpp misc_matrix_func.h
	g++ -c misc_matrix_func.cpp $(CFLAGS) $(OPTI_FLAGS)

