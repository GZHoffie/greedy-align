#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <functional>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/all.hpp>

class KPolymer {


public:
    std::vector<unsigned int> counts;
    std::vector<seqan3::dna4> nucleotides;
    unsigned int K;

public:
    KPolymer(unsigned int k, std::vector<unsigned int> c, std::vector<seqan3::dna4> n) : K(k), counts(c), nucleotides(n) {}

    std::string to_string() const {
        std::string res;
        for (size_t i = 0; i < counts.size(); ++i) {
            res += std::to_string(counts[i]) + nucleotides[i].to_char();
        }
        return res;
    }

    unsigned int getK() const {
        return K;
    }

    bool operator==(const KPolymer& other) const {
        return counts == other.counts && nucleotides == other.nucleotides;
    }
};

// Hash function for KPolymer
namespace std {
    template <>
    struct hash<KPolymer> {
        size_t operator()(const KPolymer& kPolymer) const {
            size_t hashValue = 0;
            const auto& k = kPolymer.getK();
            for (size_t i = 0; i < k; ++i) {
                for (size_t j = 0; j < kPolymer.counts[i]; ++j) {
                    hashValue |= kPolymer.nucleotides[i].to_rank();
                    hashValue <<= 2;
                }
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

    // Convert the sequence into a vector of polymers
    std::vector<unsigned int> counts;
    std::vector<seqan3::dna4> nucleotides;

    for (size_t i = 0; i < sequence.size(); ++i) {
        if (i == 0 || sequence[i] != sequence[i - 1]) {
            counts.push_back(1);
            nucleotides.push_back(sequence[i]);
        } else {
            ++counts.back();
        }
    }

    // Create KPolymer objects
    if (counts.size() < k) {
        KPolymer k_polymer(counts.size(), counts, nucleotides);
        k_polymers.push_back(k_polymer);
        return k_polymers;
    }

    for (size_t i = 0; i < counts.size() - k + 1; ++i) {
        std::vector<unsigned int> c(counts.begin() + i, counts.begin() + i + k);
        std::vector<seqan3::dna4> n(nucleotides.begin() + i, nucleotides.begin() + i + k);
        KPolymer k_polymer(k, c, n);
        k_polymers.push_back(k_polymer);

    }

    return k_polymers;
}
