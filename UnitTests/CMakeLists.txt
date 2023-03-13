include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(test_fmindex
	${CMAKE_CURRENT_SOURCE_DIR}/unit_tests.cpp
)

target_include_directories(test_fmindex
  PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/
)

target_link_libraries(test_fmindex
  GTest::gtest_main
  gmock
  FMIndex
  BWT
)

include(GoogleTest)


gtest_add_tests(TARGET test_fmindex
	TEST_SUFFIX .noArgs
	TEST_LIST   noArgsTests
)

#gtest_discover_tests(test_fmindex
	##TEST_FILTER *MixColumns*
#)