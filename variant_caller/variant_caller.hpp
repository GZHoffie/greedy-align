#ifndef GREEDY_ALIGN_VARIANT_CALLER_HPP
#define GREEDY_ALIGN_VARIANT_CALLER_HPP


#include "../utils.hpp"

#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>

class variant_caller {
protected:

public:
    variant_caller(){};

    std::string name;

    /**
     * @brief Data structure containing alignment results.
     * @param score alignment score, used to indicate mapping significance.
     * @param CIGAR the CIGAR string to indicate how the two strings are aligned together.
     */
    typedef struct {
        int score = 0;
        std::string CIGAR;
    } align_result_t;

    void read_reference(const seqan3::dna4_vector& ref);

    void process_read(const seqan3::dna4_vector& read, unsigned int position);

    std::vector<unsigned int> variant_positions();

    /**
     * @brief Virtual function that aligns two seqan3::dna4 vectors together.
     * @param s1, s2 the DNA strings to be aligned together.
     * @return align_result_t
     */
    virtual align_result_t align(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) = 0;
};

#endif