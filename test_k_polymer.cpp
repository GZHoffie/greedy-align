#include "./type/k_polymer.hpp"

using namespace seqan3::literals;
// Example usage
int main() {
    seqan3::dna4_vector sequence = "CCAATTTGGA"_dna4;
    size_t k = 3;

    auto k_polymers = to_k_polymers(sequence, k);

    for (const auto& k_polymer : k_polymers) {
        std::cout << k_polymer.to_string() << std::endl;
        std::cout << std::endl;
    }

    return 0;
}
