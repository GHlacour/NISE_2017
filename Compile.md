# This file contanins examples of how the program should be compiled on different systems
Also check the manual

## MacOS 15.0.1 with homebrew
brew install libomp
brew install imagemagick

mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER="clang" -DOpenMP_C_LIB_NAMES="libomp" -DOpenMP_CXX_LIB_NAMES="libomp" -DOpenMP_libomp_LIBRARY="/opt/homebrew/opt/libomp/lib/libomp.dylib" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
make
