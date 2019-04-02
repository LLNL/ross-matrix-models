PROJECT(MatrixCancel CXX)

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}")

LIST(APPEND mc_srcs
  matrix-cancel-model.cc)

ADD_EXECUTABLE(matrix-cancel ${mc_srcs})
TARGET_LINK_LIBRARIES(matrix-cancel ROSS)
