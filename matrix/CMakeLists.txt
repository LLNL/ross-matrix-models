PROJECT(MatrixModel CXX)

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}")

LIST(APPEND mc_srcs
  matrix-model.cc)

ADD_EXECUTABLE(matrix-model ${mc_srcs})
TARGET_LINK_LIBRARIES(matrix-model ROSS)
