PROJECT(MatrixModel CXX)

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}")

LIST(APPEND mc_srcs
  matrix-model.cc)

FIND_PACKAGE(MPI)
ADD_EXECUTABLE(matrix-model ${mc_srcs})
TARGET_LINK_LIBRARIES(matrix-model ROSS ${MPI_CXX_LIBRARIES})

ADD_TEST(matrix_model_test mpirun -np 5 ./matrix-model --synch=3 --nlp=2 --end=1000.0 --extramem=1000 --report-interval=0.1)
