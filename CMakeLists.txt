PROJECT(PointCloundNormals)
cmake_minimum_required(VERSION 3.11)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
set (CMAKE_CXX_STANDARD 11)

find_package(CGAL COMPONENTS Core Eigen3 BLAS LAPACK)
IF(CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
  ADD_DEFINITIONS("-DCGAL_EIGEN3_ENABLED   ")
  ADD_DEFINITIONS("-DWITH_CGAL ")
  SET(Deps ${Deps} ${CGAL_LIBRARIES} ${CGAL_3D_PARTY-LIBRARIES} )
  ## Making sure that CGAL got the Eigen3 flag
  ADD_DEFINITIONS("-DWITH_Eigen3 -DWITH_LAPACK ")
  message(STATUS "CGAL found.")
ENDIF(CGAL_FOUND)

SET(SRCs
  normalsPCA
  normalsJet
  normalsVCM)


FOREACH(FILE ${SRCs})
  add_executable(${FILE} ${FILE})
target_link_libraries (${FILE}  ${Deps} )
ENDFOREACH(FILE)
