# set minimum cmake version
cmake_minimum_required( VERSION 3.10.2 )

# set compiler type
set ( CMAKE_CXX_COMPILER "${CXX_COMPILER_TYPE}" )

# set project name
if ( ${ACC_TYPE} MATCHES "GPU")

   message("** GPU Accelerators")

   if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18.0)
      # useful setting to compile on m100
      set ( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -mno-float128")
   endif ()

   project( LINSOL CXX CUDA )
   enable_language(CUDA)

   if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18.0)
      if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
         set(CMAKE_CUDA_ARCHITECTURES 75)
      endif()
   endif()

else()
   message("** NO Accelerators")
   project( DRIVER CXX )
endif()

# set names
set( CORE_NAME "driver_Lagrange" )

# set install directory
set( CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} )

# set build type (DEBUG/RELEASE)
if ( ${BUILD_TYPE} MATCHES "DEBUG" )
   message("** Compiling in debug mode")
   set ( CMAKE_BUILD_TYPE Debug )
elseif ( ${BUILD_TYPE} MATCHES "RELEASE"   )
   message("** Compiling in release mode")
   set ( CMAKE_BUILD_TYPE RelWithDebInfo )
else ()
   message("+++++++++++++++++++++++++++++++++++++++++++")
   message("++ ERROR: BUILD_TYPE NOT YET IMPLEMENTED ++")
   message("+++++++++++++++++++++++++++++++++++++++++++")
endif()

# OMP
FIND_PACKAGE(OpenMP REQUIRED)

# compiler setting
set (CMAKE_CXX_FLAGS_DEBUG          "-O0 -g -std=c++11 ${OpenMP_CXX_FLAGS} -Wall")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3    -std=c++11 ${OpenMP_CXX_FLAGS} -Wall")
if ( ${ACC_TYPE} MATCHES "GPU")
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -DGPU")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DGPU")
endif()

# compiler setting for GPU
if ( ${ACC_TYPE} MATCHES "GPU")
   # old fry architecture
   set (CMAKE_CUDA_FLAGS_DEBUG          "-O0 -g -arch=sm_37 -use_fast_math -ccbin=g++ --relocatable-device-code=true -Xcompiler ${OpenMP_CXX_FLAGS} --compiler-options -Wall")
   set (CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O3    -arch=sm_37 -use_fast_math -ccbin=g++ --relocatable-device-code=true -Xcompiler ${OpenMP_CXX_FLAGS} --compiler-options -Wall")
endif()

# MPI
find_package(MPI REQUIRED COMPONENTS CXX)
include_directories(${MPI_INCLUDE_PATH})

# compiler setting for MPI
set (CMAKE_MPI_CXX_FLAGS_DEBUG          "-O0 -g -std=c++11 ${OpenMP_CXX_FLAGS} -Wall")
set (CMAKE_MPI_CXX_FLAGS_RELWITHDEBINFO "-O3    -std=c++11 ${OpenMP_CXX_FLAGS} -Wall")

# set MPI THREAD communications (FUNNELED/MULTIPLE)
if ( ${MPI_THREAD_COMM_TYPE} MATCHES "FUNNELED")
   message("** Using the MPI thread funneled mode")
elseif ( ${MPI_THREAD_COMM_TYPE} MATCHES "MULTIPLE")
   message("** Using the MPI thread multiple mode")
   message("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
   message("++ ERROR: MPI_THREAD_COMM_TYPE NOT YET IMPLEMENTED ++")
   message("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -DMultThrCom")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DMultThrCom")
else ()
   message("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
   message("++ ERROR: MPI_THREAD_COMM_TYPE NOT YET IMPLEMENTED ++")
   message("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
endif()

# find packages
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/FindPACK)
## CHRONOS
find_package( CHRONOS REQUIRED )
include_directories(${CHRONOS_INCLUDE_DIR})

# set LAPACKE or MKL
if ( ${LINEAR_ALG_LIB_TYPE} MATCHES "LAPACKE" )
   message("** Using the LAPACKE lib")
   if ( ${SPMV_KERNEL_TYPE} MATCHES "MKL_SPMV" )
      message("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
      message("++ ERROR: LINEAR_ALG_LIB_TYPE MUST BE MKL TO USE MKL_SPMV ++")
      message("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
   endif()
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -DLAPACKE")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DLAPACKE")
   find_package( USERLAPACK REQUIRED )
   include_directories(${USERLAPACK_INCLUDE_DIR})
   if ( ${Fortran_COMPILER_TYPE} MATCHES "gfortran" )
      set( LINKER_FORTRAN "-lgfortran")
   elseif ( ${Fortran_COMPILER_TYPE} MATCHES "ifort" )
      set( LINKER_FORTRAN "-lifcore")
   endif()
   set( LINKER_LIBS ${USERLAPACKE_LIB} ${USERCBLAS_LIB} ${USERLAPACK_LIB}
        ${USERBLAS_LIB} ${LINKER_FORTRAN} )
elseif ( ${LINEAR_ALG_LIB_TYPE} MATCHES "MKL" AND ${SPMV_KERNEL_TYPE} MATCHES "M3E_SPMV" )
   message("** Using the MKL lib")
   message("** Using the M3E SPMV kernel")
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -DMKL")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DMKL")
   set ( MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
   include_directories(${MKL_INCLUDE_DIR})
   set( MKL_LINKER_FLAGS "-mkl" )
   set( LINKER_LIBS ${MKL_LINKER_FLAGS} )
elseif ( ${LINEAR_ALG_LIB_TYPE} MATCHES "MKL" AND ${SPMV_KERNEL_TYPE} MATCHES "MKL_SPMV" )
   message("** Using the MKL lib")
   message("** Using the MKL SPMV kernel")
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -DMKL -DMKL_SPMV")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DMKL -DMKL_SPMV")
   set ( MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
   include_directories(${MKL_INCLUDE_DIR})
   set( MKL_LINKER_FLAGS "-mkl" )
   set( LINKER_LIBS ${MKL_LINKER_FLAGS} )
else ()
   message("++++++++++++++++++++++++++++++++++++++++++++++++++++")
   message("++ ERROR: LINEAR_ALG_LIB_TYPE OR SPMV_KERNEL_TYPE ++")
   message("++        NOT YET IMPLEMENTED                     ++")
   message("++++++++++++++++++++++++++++++++++++++++++++++++++++")
endif()

# MPI
set( LINKER_LIBS ${LINKER_LIBS} ${MPI_CXX_LIBRARIES} )

# CUB
if ( ${ACC_TYPE} MATCHES "GPU")
   message("** Using the CUB lib")
   find_package( CUB REQUIRED )
   include_directories(${CUB_INCLUDE_DIR})
endif()

# CUBLAS
if ( ${ACC_TYPE} MATCHES "GPU")
   message("** Using the CUBLAS lib")
   set( LINKER_LIBS ${LINKER_LIBS} "-lcublas" )
endif()


# METIS
message("** Using the METIS lib")
find_package( METIS REQUIRED )
include_directories(${PARMETIS_INCLUDE_DIR})
set( LINKER_LIBS ${LINKER_LIBS} ${PARMETIS_LIB} ${METIS_LIB} )

# PUGIXML
message("** Using the PUGIXML lib")
find_package( PUGIXML REQUIRED )
include_directories(${PUGIXML_INCLUDE_DIR})
set( LINKER_LIBS ${LINKER_LIBS} ${PUGIXML_LIB} )

# LEXACTIVATOR
message("** Using the LEXACTIVATOR lib")
find_package( LEXACTIVATOR REQUIRED )
include_directories(${LEXACTIVATOR_INCLUDE_DIR})
set( LINKER_LIBS ${LINKER_LIBS} ${LEXACTIVATOR_LIB} )

# JWT
message("** Using the JWT lib")
find_package( JWT REQUIRED )
include_directories(${JWT_INCLUDE_DIR})

# CURLPP
message("** Using the CURLPP lib")
find_package( CURLPP REQUIRED )
include_directories(${CURLPP_INCLUDE_DIR})
set( LINKER_LIBS ${LINKER_LIBS} ${CURLPP_LIB} )

#OpenSSL Crypto
message("** Using the OpenSSL crypto lib")
find_package( OpenSSL REQUIRED )
if (${OPENSSL_FOUND}) 
   include_directories(${OPENSSL_INCLUDE_DIR})
   set( LINKER_LIBS ${LINKER_LIBS} ${OPENSSL_CRYPTO_LIBRARIES} )
else()
   message("++++++++++++++++++++++++++++++++++++++++++++++++++++")
   message("++ ERROR: OpenSSL crypto not found                ++")
   message("++++++++++++++++++++++++++++++++++++++++++++++++++++")
endif()

# curl
message("** Using the curl lib")
find_package( CURL REQUIRED )
if (${CURL_FOUND}) 
   include_directories(${CURL_INCLUDE_DIRS})
   set( LINKER_LIBS ${LINKER_LIBS} ${CURL_LIBRARIES} )
else()
   message("++++++++++++++++++++++++++++++++++++++++++++++++++++")
   message("++ ERROR: curl not found                          ++")
   message("++++++++++++++++++++++++++++++++++++++++++++++++++++")
endif()

# set submodules
add_subdirectory( src/Matrix )
add_subdirectory( src/ReverseAugmented )
add_subdirectory( src/Utils )
add_subdirectory( src/Core )

# Dependencies
add_dependencies( Core Utils RevAugm Matrix)

