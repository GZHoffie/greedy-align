#ifndef GREEDY_ALIGN_K_MATCHING_HPP
#define GREEDY_ALIGN_K_MATCHING_HPP

#include <seqan3/search/views/kmer_hash.hpp>

#include "../utils.hpp"

template<unsigned int READ_LENGTH>
class exact_k_matching_aligner : public aligner {
private:
    std::vector<std::vector<unsigned int>> C;
    uint8_t k;

    /**
     * @brief Reset all elements in `C` to be zero.
     * 
     */
    void reset() {
        for (int i = 0; i < C.size(); i++) {
            std::fill(C[i].begin(), C[i].end(), 0);
        }
    }

    /**
     * @brief Calculate the longest common sequence between two ordered k-mer sequences. (i.e. the k-matching
     *        between the two DNA strings)
     * 
     * @param s1, s2 the two DNA strings.
     * @return the length of the longest common sequence.
     */
    unsigned int _longest_common_subsequence(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) {

        // get the k-mers of s1 and s2.
        auto k1 = s1 | seqan3::views::kmer_hash(seqan3::ungapped{k});
        auto k2 = s2 | seqan3::views::kmer_hash(seqan3::ungapped{k});

        // perform dynamic programming step
        for (int i = 1; i <= k1.size(); i++) {
            for (int j = 1; j <= k2.size(); j++) {
                if (k1[i-1] == k2[j-1]) {
                    // if the last two k-mers match
                    C[i][j] = C[i-1][j-1] + 1;
                } else {
                    C[i][j] = std::max(C[i-1][j], C[i][j-1]);
                }
            }
        }
        return C[k1.size()][k2.size()];
    }

public:
    exact_k_matching_aligner(uint8_t seed_length) {
        k = seed_length;
        for (int i = 0; i < READ_LENGTH + 1; i++) {
            std::vector<unsigned int> C_i;
            for (int j = 0; j < READ_LENGTH + 1; j++) {
                C_i.push_back(0);
            }
            C.push_back(C_i);
        }
    }

    align_result_t align(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) {
        align_result_t res;
        res.score = _longest_common_subsequence(s1, s2);
        return res;
    }
};

#endif