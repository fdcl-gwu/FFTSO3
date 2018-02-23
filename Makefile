INCLUDE_PATH= ./eigen-3.3.4
CFLAGS=$(foreach d, $(INCLUDE_PATH), -I$d) -Wall -std=c++11
VRPN_LIBS = -lvrpn -lquat -pthread
OPTI_FLAG = -O3
GTK_LIBS = -rdynamic `pkg-config --cflags gtk+-3.0`  `pkg-config --libs gtk+-3.0` 

test_FFTSO3: test_FFTSO3.cpp
	g++ -o test_FFTSO3  test_FFTSO3.cpp -lm $(CFLAGS) $(OPTI_FLAGS)

test_fdcl_nn.o: test_fdcl_nn.cpp fdcl_nn.o
	g++ -c test_fdcl_nn.cpp -lm $(CFLAGS)

fdcl_nn.o: fdcl_nn.h fdcl_nn.cpp fdcl_softmax_layer.o fdcl_mlp_layer.o fdcl_layer.o
	g++ -c fdcl_nn.cpp $(CFLAGS)

fdcl_softmax_layer.o: fdcl_softmax_layer.cpp fdcl_softmax_layer.h fdcl_layer.o
	g++ -c fdcl_softmax_layer.cpp $(CFLAGS)

fdcl_mlp_layer.o: fdcl_mlp_layer.cpp fdcl_mlp_layer.h fdcl_layer.o
	g++ -c fdcl_mlp_layer.cpp $(CFLAGS)

fdcl_layer.o: fdcl_layer.cpp fdcl_layer.h
	g++ -c fdcl_layer.cpp $(CFLAGS)
	
