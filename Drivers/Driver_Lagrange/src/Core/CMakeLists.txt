#set sources
file( GLOB FILES *.cpp )
set( CORE_SOURCES ${FILES} )
include_directories( ${PROJECT_SOURCE_DIR}/src/Utils/include)
include_directories( ${PROJECT_SOURCE_DIR}/src/Matrix/include)
include_directories( ${PROJECT_SOURCE_DIR}/src/ReverseAugmented/include)
include_directories( ${PROJECT_SOURCE_DIR}/src/Core/include)
include_directories( ${CHRONOS_INCLUDE_DIR}) 

# create executable
add_executable( Core ${CORE_SOURCES} )
set_target_properties(Core PROPERTIES OUTPUT_NAME ${CORE_NAME})
get_directory_property( RevAugm_OBJS   DIRECTORY ${PROJECT_SOURCE_DIR}/src/ReverseAugmented DEFINITION list_OBJS )
get_directory_property( Matrix_OBJS   DIRECTORY ${PROJECT_SOURCE_DIR}/src/Matrix DEFINITION list_OBJS )

target_link_libraries ( Core Utils ${RevAugm_OBJS} ${Matrix_OBJS} ${CHRONOS_LIBRARIES} ${LINKER_LIBS} )
