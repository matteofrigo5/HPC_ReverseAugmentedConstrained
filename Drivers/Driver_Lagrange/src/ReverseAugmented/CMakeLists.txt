# set all sources
file( GLOB FILES *.cpp )
set( RevAugm_SOURCES ${FILES} )


# set all headers
file( GLOB HEADERS ${PROJECT_SOURCE_DIR}/src/Preconditioner/ReverseAugmented/include/*.h )
set( RevAugm_HEADERS ${HEADERS} )


# include directories
include_directories( ${PROJECT_SOURCE_DIR}/src/Matrix/include
                     ${PROJECT_SOURCE_DIR}/src/ReverseAugmented/include )


# add obj
add_library( RevAugm STATIC OBJECT ${RevAugm_SOURCES} )

# set list of objects
set( list_OBJS $<TARGET_OBJECTS:RevAugm> )
set( list_HDRS ${RevAugm_HEADERS} )



