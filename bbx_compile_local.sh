# These are compilation commands that should work in a typical environment. 

# This step is needed only the first time after the code has been downloaded. 
# It creates the configure script suitable for your environment. 
# For subsequent recompilations this may be omitted
autoreconf -fi

# Set the C compiler, expect an MPI wrapper.
# This could be omitted if you don't have mpicc or you only want the sequential version.
export CC=mpicc

# You don't need to repeat this step if you simply recompiling the program after modifying the code. 
# Configure the system, assuming you don't have system privileges and install the binaries into your ~/bin/ directory. 
# Naturally, the ~/bin directory should exist. 
./configure CFLAGS="-O3" --prefix=$HOME
# The variant below is you want to be able to debug the executable
#./configure CFLAGS="-g -O0" --prefix=$HOME

# Make the distribution.
# This will recompile the bits that need recompiling
make

# Install the binaries, in your $HOME/bin directory,
# or whatever you specified at configure step. 
make install

