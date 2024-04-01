

#include "aligner/greedy_aligner.hpp"
#include "./utils.hpp"
#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3::literals;

int main() {



    greedy_aligner<160> ga(3, 2, 2, true);
    seqan3::dna4_vector reference = "ACGAGTCCACT"_dna4;
    seqan3::dna4_vector query =     "ACTAGAACT"_dna4;
    auto align_res = ga.align(reference, query);


    auto alignment = seqan3::alignment_from_cigar(align_res.cigar.get_cigar(), reference, 0, query);
    seqan3::debug_stream << alignment << "\n";
    
}