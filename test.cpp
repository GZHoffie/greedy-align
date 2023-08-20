#include "aligner/greedy_aligner.hpp"
#include "./utils.hpp"

int main() {

    using namespace seqan3::literals;
    greedy_aligner<160> ga(3, 2, 2, true);
    ga.align("ACTAAGATCA"_dna4,
             "ACTAACATCA"_dna4);
    seqan3::dna4_vector test("ACTAAGATCA"_dna4);
    auto vec = test | seqan3::views::to_char;
    std::string str(vec.begin(), vec.end());
    std::cout << str << "\n";
    return 0;
}
