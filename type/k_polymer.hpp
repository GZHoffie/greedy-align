#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <functional>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/dna4_vector.hpp>

class KPolymer {
public:
    using PolymerElement = std::tuple<int, char>; // (count, nucleotide)
    using PolymerVector = std::vector<PolymerElement>;

private:
    PolymerVector polymer;
    unsigned int K;

public:
    KPolymer(unsigned int k) : K(k) {}


    KPolymer(const PolymerVector& elements) : polymer(elements) {}

    const PolymerVector& getPolymer() const {
        return polymer;
    }

    bool operator==(const KPolymer& other) const {
        return polymer == other.polymer;
    }
};

// Hash function for KPolymer
namespace std {
    template <>
    struct hash<KPolymer> {
        size_t operator()(const KPolymer& kPolymer) const {
            size_t hashValue = 0;
            const auto& elements = kPolymer.getPolymer();
            for (const auto& [count, nucleotide] : elements) {
                hashValue ^= std::hash<int>()(count) ^ (std::hash<char>()(nucleotide) << 1);
            }
            return hashValue;
        }
    };
}

// Function to convert seqan3::dna4_vector to a vector of KPolymer
std::vector<KPolymer> to_k_polymers(const seqan3::dna4_vector& sequence, size_t k) {
    std::vector<KPolymer> k_polymers;

    if (k == 0 || sequence.empty()) {
        return k_polymers;
    }

    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        KPolymer::PolymerVector polymer_elements;
        size_t count = 1;
        char prev_nucleotide = seqan3::to_char(sequence[i]);

        for (size_t j = 1; j < k; ++j) {
            char current_nucleotide = seqan3::to_char(sequence[i + j]);
            if (current_nucleotide == prev_nucleotide) {
                ++count;
            } else {
                polymer_elements.emplace_back(count, prev_nucleotide);
                count = 1;
                prev_nucleotide = current_nucleotide;
            }
        }

        polymer_elements.emplace_back(count, prev_nucleotide);
        k_polymers.emplace_back(polymer_elements);
    }

    return k_polymers;
}

// Example usage
int main() {
    seqan3::dna4_vector sequence = "CCAATTTGGA"_dna4;
    size_t k = 3;

    auto k_polymers = to_k_polymers(sequence, k);

    for (const auto& k_polymer : k_polymers) {
        for (const auto& [count, nucleotide] : k_polymer.getPolymer()) {
            std::cout << count << nucleotide << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
