g++ -O2 -std=c++11 viewer.cc -framework GLUT -framework OpenGL \
    -Wno-deprecated-declarations -I$MeshFEM -I$EIGEN_INC -I$CSGFEM -I$VCGLIB_INC \
    $MeshFEM/MeshIO.o \
    -I $LIBIGL_PATH/include \
    -I $CGAL_INC -L$CGAL_LIB -lcgal -lmpfr -lgmp \
    -DHAS_TBB -I$TBB_INC -L$TBB_LIB -ltbb \
    -I /opt/local/include -L /opt/local/lib -lAntTweakBar -o viewer
