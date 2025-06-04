# SET(SRCS
# 	dgalerkin.cpp
# 	configParser.cpp
# 	Mesh.cpp
# 	utils.cpp
# 	solver.cpp
# 	eqEdit.cpp
# 	fft.cpp
# 	../include/configParser.h
# 	../include/Mesh.h
# 	../include/utils.h
# 	../include/solver.h
# 	../include/eqEdit.h
# 	../include/fft.h
# )
file(GLOB_RECURSE SRCS
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/include/*.h"
)
# 设置构建类型为 Release（开启 O3）
# set(CMAKE_BUILD_TYPE Release)

# 针对 Release 模式加优化
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -fopenmp")

# 针对 Debug 模式加调试选项（一般默认包含 -g）
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")

ADD_EXECUTABLE(dgalerkin ${SRCS})
TARGET_LINK_LIBRARIES(dgalerkin ${GMSH_LIBRARIES} ${LAPACKBLAS_LIBRARIES} ${VTK_LIBRARIES}  -fopenmp -lpthread -lwave -lfftw3) #-ltbb

vtk_module_autoinit(
	TARGETS dgalerkin
	MODULES ${VTK_LIBRARIES}
)

install(TARGETS dgalerkin
  RUNTIME DESTINATION bin
)
