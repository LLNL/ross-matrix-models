PROJECT(MatrixCancel CXX)

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}")

LIST(APPEND mc_srcs
  matrix-cancel-model.cc)

FIND_PACKAGE(MPI)
ADD_EXECUTABLE(matrix-cancel ${mc_srcs})
TARGET_LINK_LIBRARIES(matrix-cancel ROSS ${MPI_CXX_LIBRARIES})

ADD_TEST(matrix_cancel_serial ./matrix-cancel --nlp=10240 --end=100 --synch=1 --extramem=1000 --report-interval=0.1 --batch=2 --gvt-interval=1024)
FOREACH(cpn 3 4 5 8 13 16)
  MATH(EXPR nlp "10240 / ${cpn}")
  ADD_TEST(matrix_cancel_${cpn} mpirun -np ${cpn} ./matrix-cancel --nlp=${nlp} --end=100 --synch=3 --extramem=1000 --report-interval=0.1 --batch=2 --gvt-interval=1024)
ENDFOREACH(cpn)
