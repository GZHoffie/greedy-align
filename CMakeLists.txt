cmake_minimum_required(VERSION 3.4)
project(greedy-aligner CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")


include(FetchContent)

# seqan3 for the basic skeleton for the mapper
FetchContent_Declare(seqan3
    GIT_REPOSITORY https://github.com/seqan/seqan3.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(seqan3)

# Sharg for argument parsing
FetchContent_Declare(sharg
    GIT_REPOSITORY https://github.com/seqan/sharg-parser.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(sharg)

#FetchContent_Declare(parasail
#  GIT_REPOSITORY https://github.com/jeffdaily/parasail.git
#  GIT_TAG        master
#)
#FetchContent_MakeAvailable(parasail)

# set common header files
set(ALIGNER_FILES
    aligner/aligner.hpp
    aligner/sw_aligner.hpp
    aligner/greedy_aligner.hpp
    aligner/k_matching.hpp)


set(UTIL_FILES
    utils.hpp
    tools/generate_dataset.hpp)

set(SHARED_FILES
    ${UTIL_FILES}
    ${ALIGNER_FILES})

# build dataset generating tools
add_executable(greedy_aligner main.cpp ${SHARED_FILES})
target_link_libraries(greedy_aligner PUBLIC seqan3::seqan3 sharg::sharg)


add_executable(greedy_aligner_test test.cpp ${SHARED_FILES})
target_link_libraries(greedy_aligner_test PUBLIC seqan3::seqan3)
