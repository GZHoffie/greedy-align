#include "aligner/greedy_aligner.hpp"
#include "./utils.hpp"

int main() {

    using namespace seqan3::literals;
    greedy_aligner<10> ga(3, 2);
    ga.align("ACGTGCA"_dna4, "ACCTGCA"_dna4);
    return 0;
}