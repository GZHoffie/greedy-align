#ifndef GREEDY_ALIGN_PARASAIL_ALIGNER
#define GREEDY_ALIGN_PARASAIL_ALIGNER

#include "./aligner.hpp"

#include "parasail.h"

class parasail_aligner : public aligner {
private:
    // SIMD Needleman-Wunsch objects
    parasail_matrix_t* penalty_matrix;

    /**
     * @brief Convert a seqan3::dna4_vector object to a string.
     * 
     * @param s the vector of seqan3::dna4 objects.
     */
    const char* _to_string(const seqan3::dna4_vector& s) {
        auto vec = s | seqan3::views::to_char;
        std::string str(vec.begin(), vec.end());
        return str.c_str();
    }

public:
    parasail_aligner() {
        penalty_matrix = parasail_matrix_create("ACGT", 0, -1);
    }

    align_result_t align(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) {
        // initialize result
        align_result_t res;
        
        // convert both to char*
        // TODO: find more efficient ways of conversion.
        parasail_result_t* result;
        parasail_cigar_t* cigar_result;

        // perform alignment
        const char* c1 = _to_string(s1);
        const char* c2 = _to_string(s2);
        result =  parasail_nw_trace_striped_avx2_256_16(c1, s1.size(), c2, s2.size(), 1, 1, penalty_matrix);
        cigar_result = parasail_result_get_cigar(result, c1, s1.size(), c2, s2.size(), penalty_matrix);
        res.CIGAR = parasail_cigar_decode(cigar_result);
        res.score = result->score;
        
        // free allocated memory
        parasail_result_free(result);
        parasail_cigar_free(cigar_result);

        return res;
    }
};

#endif