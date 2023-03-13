add_library(FMIndex)

target_sources(FMIndex
    PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/BitVector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FMIndex.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/WaveletTree.cpp
)

target_include_directories(FMIndex
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/
	${PROJECT_BINARY_DIR}
) 

target_link_libraries(FMIndex
  BWT
)
