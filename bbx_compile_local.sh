# These are compilation commands that should work in a typical environment. 

# This step is needed only when you have made substantial changes to BeatBox, including
# creating new source files etc. 
# It re-creates the configure script suitable for your environment. 
# For compilation of the downloaded BeatBox distribution with a ready configure script,
# this is not needed.
# Uncomment the next line only if you have a working version of autotools installed.
# autoreconf -fi

# Set the C compiler, expect an MPI wrapper.
# This could be omitted if you don't have mpicc or you only want the sequential version.
# In that case, comment out the next line.  
export CC=mpicc

# You don't need to repeat the next step if you simply recompiling the program after modifying the code. 
# Configure the system, assuming you don't have system privileges and install the binaries
# into your ~/bin/ directory. 
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

