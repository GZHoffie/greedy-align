cmake_minimum_required(VERSION 3.4)
project(greedy-aligner CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")


include(FetchContent)

# seqan3 for the basic skeleton for the mapper
FetchContent_Declare(seqan3
    GIT_REPOSITORY https://github.com/seqan/seqan3.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(seqan3)

# Sharg for argument parsing
FetchContent_Declare(sharg
    GIT_REPOSITORY https://github.com/seqan/sharg-parser.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(sharg)

FetchContent_Declare(parasail
  GIT_REPOSITORY https://github.com/jeffdaily/parasail.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(parasail)

# set common header files
set(ALIGNER_FILES
    aligner/aligner.hpp
    aligner/greedy_aligner.hpp
    ensembler/ensembler.hpp
    ensembler/pog_ensembler.hpp)

set(VARIANT_CALLER_FILES
    variant_caller/variant_caller.hpp
    variant_caller/greedy_variant_caller.hpp)

set(BENCHMARK_FILES
    aligner/seqan_aligner.hpp
    aligner/parasail_aligner.hpp
    aligner/k_matching.hpp)

set(UTIL_FILES
    utils.hpp
    tools/pairwise_dataset.hpp
    tools/msa_dataset.hpp
    benchmark/benchmark.hpp)

set(SHARED_FILES
    ${UTIL_FILES}
    ${ALIGNER_FILES})

# build dataset generating tools
add_executable(greedy_aligner main.cpp ${SHARED_FILES})
target_link_libraries(greedy_aligner PUBLIC seqan3::seqan3 sharg::sharg parasail)


add_executable(greedy_aligner_test test.cpp ${SHARED_FILES})
target_link_libraries(greedy_aligner_test PUBLIC seqan3::seqan3)

add_executable(aligner_benchmark benchmark.cpp ${BENCHMARK_FILES} ${SHARED_FILES})
target_link_libraries(aligner_benchmark PUBLIC seqan3::seqan3 sharg::sharg parasail)