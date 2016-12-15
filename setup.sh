# setup.sh
# Make the "myext" Python Module ("myext.so")
export CC="gcc"   
export CXX="g++"   
export CFLAGS="-I.includes/ "   
export LDFLAGS="-L."   
    python setup.py build_ext --inplace
