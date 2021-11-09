# set all sources
file( GLOB FILES *.cpp )
set( MATRIX_SOURCES ${FILES} )


# set all headers
file( GLOB HEADERS ${PROJECT_SOURCE_DIR}/src/Matrix/include/*.h )
set( MATRIX_HEADERS ${HEADERS} )


# include directories
include_directories( ${PROJECT_SOURCE_DIR}/src/Matrix/include)

# add obj
add_library( Matrix STATIC OBJECT ${MATRIX_SOURCES} )


# set list of objects
set( list_OBJS $<TARGET_OBJECTS:Matrix> )
set( list_HDRS ${MATRIX_HEADERS} )

