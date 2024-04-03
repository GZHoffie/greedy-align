

#include "aligner/greedy_aligner.hpp"
#include "./utils.hpp"
#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3::literals;

int main() {
    greedy_aligner<160> ga(4, 7, 3, true);
    ga.align("ACCATAATGCGTGGGGCCGACTCCGGAATGCGGTCTCCATGCGCGTTTCTCCAACCTAAGGTAGCCTGTAGTTCATTGACCTCTGATGGCGCTTATGAAACCGGGAA"_dna4,
             "ACCATAATGCGTGGGGCCGACCTCGGAAATGCGGTCTCCATGCGCGTTTCCTCCAACCTAAGGTAGCCTTAGGAAGTTCATTGGACTCTGATGGCGCTTATAGAAACCGGGAA"_dna4);
    seqan3::dna4_vector test("ACTAAGATCA"_dna4);
    auto vec = test | seqan3::views::to_char;
    std::string str(vec.begin(), vec.end());
    std::cout << str << "\n";
    return 0;
}
