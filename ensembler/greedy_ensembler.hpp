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

    template <typename KEY_TYPE, typename VALUE_TYPE>
    std::unordered_map<KEY_TYPE, VALUE_TYPE> top_w_count_kmers(const std::unordered_map<KEY_TYPE, VALUE_TYPE>& counts, unsigned int w) {
        if (w >= counts.size()) {
            return counts;
        }

        // find the kmers with the highest values in the `counts`
        std::vector<std::pair<KEY_TYPE, VALUE_TYPE>> count_vector(counts.begin(), counts.end());

        std::nth_element(count_vector.begin(), count_vector.begin() + w, count_vector.end(), 
            [](const std::pair<KEY_TYPE, VALUE_TYPE>& a, const std::pair<KEY_TYPE, VALUE_TYPE>& b) {
                return a.second > b.second;
            });

        std::unordered_map<KEY_TYPE, VALUE_TYPE> res;
        for (unsigned int i = 0; i < w; i++) {
            res[count_vector[i].first] = count_vector[i].second;
        }
        return res;   
    }

    template <typename KEY_TYPE, typename VALUE_TYPE>
    std::vector<std::pair<KEY_TYPE, VALUE_TYPE>> top_w_count_kmers(const std::vector<std::pair<KEY_TYPE, VALUE_TYPE>>& counts, unsigned int w, bool ascending = false) {
        if (w >= counts.size()) {
            return counts;
        }

        // find the kmers with the highest values in the `counts`
        std::vector<std::pair<KEY_TYPE, VALUE_TYPE>> count_vector(counts);

        if (ascending) {
            std::nth_element(count_vector.begin(), count_vector.begin() + w, count_vector.end(), 
                [](const std::pair<KEY_TYPE, VALUE_TYPE>& a, const std::pair<KEY_TYPE, VALUE_TYPE>& b) {
                    return a.second < b.second;
                });
        } else {
            std::nth_element(count_vector.begin(), count_vector.begin() + w, count_vector.end(), 
                [](const std::pair<KEY_TYPE, VALUE_TYPE>& a, const std::pair<KEY_TYPE, VALUE_TYPE>& b) {
                    return a.second > b.second;
                });
        }

        std::vector<std::pair<KEY_TYPE, VALUE_TYPE>> res(count_vector.begin(), count_vector.begin() + w);
        return res;   
    }


    seqan3::dna4_vector kmer_vector_to_sequence(const std::vector<unsigned int>& kmer_vector, uint8_t k) {
        seqan3::dna4_vector res;
        if (kmer_vector.size() > 0) {
            res = hash_to_kmer(kmer_vector[0], k);
            for (unsigned int i = 1; i < kmer_vector.size(); i++) {
                res.push_back(seqan3::assign_rank_to(kmer_vector[i] & 3, seqan3::dna4{}));
            }
        }
        return res;
    }

    seqan3::dna4_vector find_consensus_greedy(const std::vector<seqan3::dna4_vector>& sv, uint8_t k, unsigned int expected_length) {
        // find the initial k-mer 
        auto [initial_counts, last_counts] = initial_and_last_kmer_counts(sv, k);
        auto all_counts = all_kmer_counts(sv, k);


        // greedily find the next k-mer.
        auto current_kmer = max_count_kmer(initial_counts);
        auto target_kmer = max_count_kmer(last_counts);

        auto res = hash_to_kmer(current_kmer, k);
        
        for (unsigned int l = 0; l < expected_length - k; l++) {
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
                if (DEBUG) {
                    seqan3::debug_stream << res << "\n";
                }
                //seqan3::debug_stream << "Current: " << hash_to_kmer(current_kmer, k) << "\n";
            }
        }

        return res;
    }

    seqan3::dna4_vector find_consensus_beam(const std::vector<seqan3::dna4_vector>& sv, uint8_t k, unsigned int width, unsigned int expected_length) {
        // find the initial k-mer 
        auto [initial_counts, last_counts] = initial_and_last_kmer_counts(sv, k);
        auto all_counts = all_kmer_counts(sv, k);

        auto number_of_sequences = sv.size();

        // greedily find the next k-mer.
        auto current_kmer = max_count_kmer(initial_counts);
        auto target_kmer = max_count_kmer(last_counts);

        auto top_w_init_kmers = top_w_count_kmers<unsigned int, unsigned int>(initial_counts, width);

        // initialize variables that store the path to the k-mer
        std::vector<std::pair<std::vector<unsigned int>, unsigned int>> current_frontier;

        for (auto const [kmer, weight] : top_w_init_kmers) {
            current_frontier.push_back(std::make_pair(std::vector<unsigned int>{kmer}, weight));
        }



        // Perform greedy beam search
        for (unsigned int l = 0; l < expected_length - k; l++) {
            std::vector<std::pair<std::vector<unsigned int>, unsigned int>> new_frontier;

            for (auto [path, weight] : current_frontier) {
                auto current_kmer = path[path.size() - 1];

                // find the next k-mer with the most votes
                for (unsigned int i = 0; i < 4; i++) {
                    unsigned int next_kmer = (current_kmer << 2) & ((1 << (2*k)) - 1) | i;
                    if (all_counts.contains(next_kmer)) {
                        // update the path
                        std::vector<unsigned int> current_path(path);
                        current_path.push_back(next_kmer);
                        new_frontier.push_back(std::make_pair(current_path, weight + all_counts[next_kmer]));
                    }
                }

            }

            if (new_frontier.empty()) {
                // no more candidates
                break;
            }

            // find the top w candidates
            current_frontier = top_w_count_kmers<std::vector<unsigned int>, unsigned int>(new_frontier, width);

            if (DEBUG) {
                seqan3::debug_stream << "Frontier: ";
                for (auto [path, weight] : current_frontier) {
                    seqan3::debug_stream << kmer_vector_to_sequence(path, k) << " " << weight << "\n";
                }
                seqan3::debug_stream << "\n";
            }
        }

        // the best path is simply the longest path with the highest weight
        unsigned int max_weight = 0;
        unsigned int max_length = 0;
        std::vector<unsigned int> max_path;
        for (auto [path, weight] : current_frontier) {
            //seqan3::debug_stream << "Path: " << kmer_vector_to_sequence(path, k) << " " << weight << "\n";
            //seqan3::debug_stream << kmer_vector_to_sequence(max_path, k) << ", " <<  max_weight << "\n";
            if (path.size() > max_length || (path.size() == max_length && weight > max_weight)) {
                max_weight = weight;
                max_length = path.size();
                max_path = path;
            }
        }
        //seqan3::debug_stream << "Max kmer" << hash_to_kmer(max_kmer, k) << " " << max_weight << "\n";

        // Turn the path into a sequence
        seqan3::dna4_vector res = kmer_vector_to_sequence(max_path, k);
        //seqan3::debug_stream << "RES:" << res << "\n";


        return res;
    }

    seqan3::dna4_vector find_consensus_bfs(const std::vector<seqan3::dna4_vector>& sv, uint8_t k, unsigned int expected_length, float agreement_threshold = 0.1) {
        // find the initial k-mer 
        auto [initial_counts, last_counts] = initial_and_last_kmer_counts(sv, k);
        auto all_counts = all_kmer_counts(sv, k);

        seqan3::debug_stream << "Number of k-mers: " << all_counts.size() << "\n";

        auto number_of_sequences = sv.size();
        unsigned int agreement_threshold_count = (unsigned int)(number_of_sequences * agreement_threshold);

        // greedily find the next k-mer.
        auto current_kmer = max_count_kmer(initial_counts);
        auto target_kmer = max_count_kmer(last_counts);

        std::vector<std::tuple<unsigned int, std::vector<unsigned int>, unsigned int>> frontier;
        frontier.push_back(std::make_tuple(current_kmer, std::vector<unsigned int>{current_kmer}, 0));

        // Perform BFS
        for (unsigned int l = 0; l < expected_length - k; l++) {
            std::vector<std::tuple<unsigned int, std::vector<unsigned int>, unsigned int>> new_frontier;

            for (auto [current_kmer, path, weight] : frontier) {
                // find the next
                for (unsigned int i = 0; i < 4; i++) {
                    unsigned int next_kmer = (current_kmer << 2) & ((1 << (2*k)) - 1) | i;
                    if (all_counts.contains(next_kmer) && all_counts[next_kmer] >= agreement_threshold_count) {
                        // update the path
                        std::vector<unsigned int> current_path(path);
                        current_path.push_back(next_kmer);
                        new_frontier.push_back(std::make_tuple(next_kmer, current_path, weight + all_counts[next_kmer]));
                    }
                } 
            }

            seqan3::debug_stream << "Frontier size: " << new_frontier.size() << "\n";

            if (new_frontier.empty()) {
                // no more candidates
                break;
            }

            frontier = new_frontier;
        }

        // the best path is the path with the highest weight
        unsigned int max_weight = 0;
        std::vector<unsigned int> max_path;
        for (auto [current_kmer, path, weight] : frontier) {
            if (weight > max_weight) {
                max_weight = weight;
                max_path = path;
            }
        }

        // Turn the path into a sequence
        seqan3::dna4_vector res = kmer_vector_to_sequence(max_path, k);

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
        return find_consensus_beam(sv, k, 10, 110);
        //return find_consensus_bfs(sv, k, 110);
        //return find_consensus_greedy(sv, k, 150);
    }
};

#endif