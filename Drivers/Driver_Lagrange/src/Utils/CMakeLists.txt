# Sources
file( GLOB FILES *.cpp )

# All sources
set( UTILS_SOURCES ${FILES} )

# Includes
include_directories( ${PROJECT_SOURCE_DIR}/src/Utils/include )
include_directories( ${PROJECT_SOURCE_DIR}/src/ReverseAugmented/include )
include_directories( ${CHRONOS_INCLUDE_DIR})

# Targets
add_library( Utils ${UTILS_SOURCES} )

