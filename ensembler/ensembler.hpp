#ifndef GREEDY_ALIGN_ENSEMBLER_HPP
#define GREEDY_ALIGN_ENSEMBLER_HPP


#include "../utils.hpp"

#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>

class ensembler {
protected:

public:
    ensembler(){};

    std::string name;

    /**
     * @brief Virtual function that finds the consensus of a vector of similar DNA sequences.
     * @param sv The vector storing all the DNA sequences.
     * @return The consensus of all the DNA strings. 
     */
    virtual seqan3::dna4_vector ensemble(const std::vector<seqan3::dna4_vector>& sv) = 0;
};

#endif