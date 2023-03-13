add_library(BWT)

target_sources(BWT
    PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/BWT.c
)

target_include_directories(BWT
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/
	${PROJECT_BINARY_DIR}
) 