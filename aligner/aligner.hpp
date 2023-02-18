#ifndef GREEDY_ALIGN_ALIGNER_HPP
#define GREEDY_ALIGN_ALIGNER_HPP


#include "../utils.hpp"

#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>

class aligner {
protected:

public:
    aligner(){};

    /**
     * @brief Data structure containing alignment results.
     * @param score alignment score, used to indicate mapping significance.
     * @param CIGAR the CIGAR string to indicate how the two strings are aligned together.
     */
    typedef struct {
        int score;
        std::string CIGAR;
    } align_result_t; 

    /**
     * @brief Virtual function that aligns two seqan3::dna4 vectors together.
     * @param s1, s2 the DNA strings to be aligned together.
     * @return align_result_t
     */
    virtual align_result_t align(seqan3::dna4_vector s1, seqan3::dna4_vector s2) = 0;
};

#endif