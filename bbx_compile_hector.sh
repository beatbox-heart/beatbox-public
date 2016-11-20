# These are compilation commands that work on Hector.
# Hector is the UK National HPC service. More information
# is available from http://www.hector.ac.uk.

# This step is needed only the first time after the code has been downloaded. 
# It creates the configure script suitable for your environment. 
# For subsequent recompilations this may be omitted
autoreconf -fi

# Load the PGI compiler module: this is what makes Hector particularly special
module load pgi

# Tell the configure command what compiler we need to use.
# This is also different from typical local compile procedure. 
export CC=cc

# Need to explicitly tell the Cray system what X libraries to use
# This is also different from typical local compile procedure. 
export LIBS="-lm -lX11 -lxcb -lxcb-xlib -ldl -lXau"

# Run configure - note that it will install in $HOME/bin
# Note that on HECToR to run a program on the back-end system
# you need to be able to access the executable from your /work
# file system.
./configure --prefix=$HOME

# Now make the code
make

# Install the code to the directories specified by the prefix above.
make install

