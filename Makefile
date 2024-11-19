# Compiler and Flags
CC = gcc
CFLAGS = -O2 -I include \
			-I /usr/local/MATLAB/R2024b/extern/include
LDFLAGS = -L /usr/local/MATLAB/R2024b/bin/glnxa64 \
			-Wl,-rpath,/usr/local/MATLAB/R2024b/bin/glnxa64 \
			-lm -lpthread -lopenblas -lmat -lmx -fopenmp

# Binary 1: v0
SRC_v0 = src/mat_read_write.c src/knn_search.c src/v0.c
OBJ_v0 = $(SRC_v0:src/%.c=build/%.o)
BIN_v0 = build/v0

# Binary 2: v1
SRC_v1 = src/mat_read_write.c src/knn_search.c src/v1.c
OBJ_v1 = $(SRC_v1:src/%.c=build/%.o)
BIN_v1 = build/v1

# Binary 2: v1_openmp
SRC_v1_openmp = src/mat_read_write.c src/knn_search.c src/v1_openmp.c
OBJ_v1_openmp = $(SRC_v1_openmp:src/%.c=build/%.o)
BIN_v1_openmp = build/v1_openmp

# Default Target: Build all binaries
all: $(BIN_v0) $(BIN_v1) $(BIN_v1_openmp)

# Rule to Compile Object Files
build/%.o: src/%.c
	mkdir -p datasets
	mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to Link Binary 1
$(BIN_v0): $(OBJ_v0)
	$(CC) $(OBJ_v0) -o $@ $(LDFLAGS)

# Rule to Link Binary 2
$(BIN_v1): $(OBJ_v1)
	$(CC) $(OBJ_v1) -o $@ $(LDFLAGS)

#Rule to Link Binary 3
$(BIN_v1_openmp): $(OBJ_v1_openmp)
	$(CC) $(OBJ_v1_openmp) -o $@ $(LDFLAGS)
	rm -rf build/*.o

# Clean Target
clean:
	rm -rf build
