#ifndef GREEDY_ALIGN_GREEDY_ENSEMBLER_HPP
#define GREEDY_ALIGN_GREEDY_ENSEMBLER_HPP

#include "ensembler.hpp"


class greedy_ensembler : public ensembler {
private:
    uint8_t K_LOWER_BOUND, K_UPPER_BOUND;
    bool DEBUG;
    /**
     * Output the smallest k such that all the k-mers appear at most once 
     * in each of the input sequences.
     */
    uint8_t choose_k(const std::vector<seqan3::dna4_vector>& sv) {
        for (uint8_t res = K_LOWER_BOUND; res <= K_UPPER_BOUND; res++) {
            bool flag = true;
            for (auto s : sv) {
                auto kmers = s | seqan3::views::kmer_hash(seqan3::ungapped{(uint8_t)res});
                std::unordered_set<unsigned int> kmer_set(kmers.begin(), kmers.end());
                if (kmer_set.size() != kmers.size()) {
                    flag = false;
                    break;
                }
            }
            // return k to be one higher so that all (k-1)-mers are unique.
            if (flag) return res + 1;
        }
        return K_UPPER_BOUND;
    }

    std::unordered_map<unsigned int, unsigned int> all_kmer_counts(const std::vector<seqan3::dna4_vector>& sv, uint8_t k) {
        std::unordered_map<unsigned int, unsigned int> res;
        for (auto s : sv) {
            auto kmers = s | seqan3::views::kmer_hash(seqan3::ungapped{k});
            std::unordered_set<unsigned int> kmer_set(kmers.begin(), kmers.end());
            for (auto kmer : kmer_set) {
                if (res.contains(kmer)) {
                    res[kmer]++;
                } else {
                    res[kmer] = 1;
                }
            }
        }
        return res;
    }

    std::pair<std::unordered_map<unsigned int, unsigned int>,
              std::unordered_map<unsigned int, unsigned int>> 
    initial_and_last_kmer_counts(const std::vector<seqan3::dna4_vector>& sv, uint8_t k) {
        std::unordered_map<unsigned int, unsigned int> initial_kmer_count;
        std::unordered_map<unsigned int, unsigned int> last_kmer_count;

        // only count the first and the last k-mer in each sequence
        for (auto s : sv) {
            auto kmer = s | seqan3::views::kmer_hash(seqan3::ungapped{k});
            std::vector<unsigned int> kmer_vector(kmer.begin(), kmer.end());
            if (kmer_vector.size() > 0) {
                initial_kmer_count[kmer_vector[0]]++;
                last_kmer_count[kmer_vector[kmer_vector.size()-1]]++;
            }
        }
        return std::make_pair(initial_kmer_count, last_kmer_count);
    }

    unsigned int max_count_kmer(const std::unordered_map<unsigned int, unsigned int>& counts) {
        unsigned int max_count = 0;
        unsigned int max_kmer = 0;
        for (auto& [kmer, count] : counts) {
            if (count > max_count) {
                max_count = count;
                max_kmer = kmer;
            }
        }
        return max_kmer;
    }

    seqan3::dna4_vector find_consensus(const std::vector<seqan3::dna4_vector>& sv, uint8_t k) {
        // find the initial k-mer 
        auto [initial_counts, last_counts] = initial_and_last_kmer_counts(sv, k);
        auto all_counts = all_kmer_counts(sv, k);


        // greedily find the next k-mer.
        auto current_kmer = max_count_kmer(initial_counts);
        auto target_kmer = max_count_kmer(last_counts);

        auto res = hash_to_kmer(current_kmer, k);
        
        while (true) {
            unsigned int max_count = 0;
            unsigned int max_kmer = 0;

            // find the next k-mer with the most votes
            for (unsigned int i = 0; i < 4; i++) {
                //seqan3::debug_stream << "Current: " << current_kmer << " " << hash_to_kmer(current_kmer, k) << "\n";
                unsigned int next_kmer = (current_kmer << 2) & ((1 << (2*k)) - 1) | i;
                //seqan3::debug_stream << hash_to_kmer(next_kmer, k) << " " << all_counts[next_kmer] << "\n";
                if (all_counts.contains(next_kmer) && all_counts[next_kmer] > max_count) {
                    if (DEBUG) {
                        seqan3::debug_stream << "Next: " << hash_to_kmer(next_kmer, k) << " " << all_counts[next_kmer] << "\n";
                    }
                    max_count = all_counts[next_kmer];
                    max_kmer = next_kmer;
                }
            }
            if (DEBUG) {
                seqan3::debug_stream << "Max: " << hash_to_kmer(max_kmer, k) << " " << max_count << "\n";
            }
            if (max_count == 0) return res;
            else {
                res.push_back(seqan3::assign_rank_to(max_kmer & 3, seqan3::dna4{}));
                current_kmer = max_kmer;

                if (current_kmer == target_kmer) {
                    return res;
                }
                //seqan3::debug_stream << "Current: " << hash_to_kmer(current_kmer, k) << "\n";
            }
        }

        return res;
    }


public:
    greedy_ensembler(uint8_t lower_bound, uint8_t upper_bound, bool debug = false) {
        K_LOWER_BOUND = lower_bound;
        K_UPPER_BOUND = upper_bound;
        DEBUG = debug;
    }

    seqan3::dna4_vector ensemble(const std::vector<seqan3::dna4_vector>& sv) override {
        uint8_t k = choose_k(sv);
        if (DEBUG) {
            seqan3::debug_stream << "Chose k: " << (int)k << "\n";
        }
        //seqan3::debug_stream << "Chose k: " << (int)k << "\n";
        return find_consensus(sv, k);
    }
};

#endif