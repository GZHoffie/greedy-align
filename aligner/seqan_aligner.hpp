#ifndef GREEDY_ALIGN_SEQAN_ALIGNER
#define GREEDY_ALIGN_SEQAN_ALIGNER

#include "./aligner.hpp"


#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>

class seqan_aligner : public aligner {
private:


public:
    seqan_aligner() {
        name = "Seqan3 Aligner";
    }


    align_result_t align(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) {
        // initialize result
        align_result_t res;

        auto align_config =
            seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                             seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                             seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                             seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
            | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{}
            | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};
        
        // perform pairwise alignment
        auto results = seqan3::align_pairwise(std::tie(s1, s2), align_config);
        auto & align_result = *results.begin();
        res.score = align_result.score();
        CIGAR cigar(seqan3::cigar_from_alignment(align_result.alignment()));
        res.CIGAR = cigar.to_string();

        return res;
    }
};

#endif