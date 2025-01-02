#ifndef GREEDY_ALIGN_GREEDY_ENSEMBLER_HPP
#define GREEDY_ALIGN_GREEDY_ENSEMBLER_HPP

#include <cmath>

#include "ensembler.hpp"
#include "../type/k_polymer.hpp"

template <typename KMER_TYPE>
class greedy_ensembler : public ensembler {
public:
    using kmer_counter_t = std::unordered_map<KMER_TYPE, unsigned int>;
    using kmer_set_t = std::unordered_set<KMER_TYPE>;
    using kmer_vector_t = std::vector<KMER_TYPE>;
    using neighbor_map_t = std::unordered_map<KMER_TYPE, std::vector<KMER_TYPE>>;

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
                kmer_vector_t kmers;
                if (std::is_same<KMER_TYPE, KPolymer>::value) {
                    auto kmer_range = to_k_polymers(s, res);
                    //kmers = kmer_vector_t(kmer_range.begin(), kmer_range.end());
                } else {
                    auto kmer_range = s | seqan3::views::kmer_hash(seqan3::ungapped{(uint8_t)res});
                    kmers = kmer_vector_t(kmer_range.begin(), kmer_range.end());
                }
                //auto kmers = s | seqan3::views::kmer_hash(seqan3::ungapped{(uint8_t)res});
                kmer_set_t kmer_set(kmers.begin(), kmers.end());
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

    kmer_counter_t all_kmer_counts(const std::vector<seqan3::dna4_vector>& sv, uint8_t k) {
        kmer_counter_t res;
        for (auto s : sv) {
            auto kmer_range = s | seqan3::views::kmer_hash(seqan3::ungapped{k});
            auto kmers = kmer_set_t(kmer_range.begin(), kmer_range.end());
            for (auto kmer : kmers) {
                if (res.contains(kmer)) {
                    res[kmer]++;
                } else {
                    res[kmer] = 1;
                }
            }
        }
        return res;
    }

    neighbor_map_t all_neighbors(const kmer_counter_t& counts) {
        neighbor_map_t res;
        for (auto& [kmer, count] : counts) {
            seqan3::dna4_vector k_1_mer_nucleotides(kmer.nucleotides.begin() + 1, kmer.nucleotides.end());
            std::vector<unsigned int> k_1_mer_counts(kmer.counts.begin() + 1, kmer.counts.end());
            KMER_TYPE previous_k_1_mer(kmer.getK() - 1, k_1_mer_counts, k_1_mer_nucleotides);
            if (res.contains(previous_k_1_mer)) {
                res[previous_k_1_mer].push_back(kmer);
            } else {
                res[previous_k_1_mer] = std::vector<KMER_TYPE>{kmer};
            }
        }
        return res;
    }

    std::pair<kmer_counter_t, kmer_counter_t> 
    initial_and_last_kmer_counts(const std::vector<seqan3::dna4_vector>& sv, uint8_t k) {
        kmer_counter_t initial_kmer_count;
        kmer_counter_t last_kmer_count;

        // only count the first and the last k-mer in each sequence
        for (auto s : sv) {
            //auto kmer = s | seqan3::views::kmer_hash(seqan3::ungapped{k});
            kmer_vector_t kmers;
            if (std::is_same<KMER_TYPE, KPolymer>::value) {
                auto kmer_range = to_k_polymers(s, k);
                //kmers = kmer_vector_t(kmer_range.begin(), kmer_range.end());
            } else {
                auto kmer_range = s | seqan3::views::kmer_hash(seqan3::ungapped{k});
                kmers = kmer_vector_t(kmer_range.begin(), kmer_range.end());
            }
            std::vector<KMER_TYPE> kmer_vector(kmers.begin(), kmers.end());
            if (kmer_vector.size() > 0) {
                initial_kmer_count[kmer_vector[0]]++;
                last_kmer_count[kmer_vector[kmer_vector.size()-1]]++;
            }
        }
        return std::make_pair(initial_kmer_count, last_kmer_count);
    }

    KMER_TYPE max_count_kmer(const kmer_counter_t& counts) {
        unsigned int max_count = 0;
        KMER_TYPE max_kmer;
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


    seqan3::dna4_vector kmer_vector_to_sequence(const std::vector<KMER_TYPE>& kmer_vector, uint8_t k) {
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
            KMER_TYPE max_kmer = 0;

            // find the next k-mer with the most votes
            for (unsigned int i = 0; i < 4; i++) {
                //seqan3::debug_stream << "Current: " << current_kmer << " " << hash_to_kmer(current_kmer, k) << "\n";
                KMER_TYPE next_kmer = (current_kmer << 2) & ((1 << (2*k)) - 1) | i;
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

    std::vector<std::pair<std::vector<unsigned int>, float>> 
    beam_search_step(const std::vector<std::pair<std::vector<unsigned int>, float>>& current_frontier, const std::unordered_map<unsigned int, unsigned int>& all_counts, unsigned int width, bool forward, uint8_t k) {
        std::vector<std::pair<std::vector<unsigned int>, float>> new_frontier;

        for (auto [path, weight] : current_frontier) {
            auto current_kmer = path[path.size() - 1];

            // find the next k-mer with the most votes
            for (unsigned int i = 0; i < 4; i++) {
                unsigned int next_kmer;
                if (forward) {
                    next_kmer = (current_kmer << 2) & ((1 << (2*k)) - 1) | i;
                } else {
                    next_kmer = (current_kmer >> 2) | (i << (2*k - 2));
                }
                if (all_counts.contains(next_kmer)) {
                    // update the path
                    std::vector<unsigned int> current_path(path);
                    current_path.push_back(next_kmer);
                    new_frontier.push_back(std::make_pair(current_path, weight + std::log2((float)all_counts.at(next_kmer))));
                }
            }
        }

        if (new_frontier.empty()) {
            // no more candidates
            return new_frontier;
        }

        // find the top w candidates
        return top_w_count_kmers<std::vector<unsigned int>, float>(new_frontier, width);
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
        std::vector<std::pair<std::vector<unsigned int>, float>> current_frontier;

        //for (auto const [kmer, weight] : top_w_init_kmers) {
        //    current_frontier.push_back(std::make_pair(std::vector<unsigned int>{kmer}, std::log2((float)weight)));
        //}

        current_frontier.push_back(std::make_pair(std::vector<unsigned int>{current_kmer}, std::log2((float)initial_counts[current_kmer])));



        // Perform greedy beam search
        for (unsigned int l = 0; l < expected_length - k; l++) {
            auto new_frontier = beam_search_step(current_frontier, all_counts, width, true, k);

            if (new_frontier.empty()) {
                // no more candidates
                break;
            }

            // find the top w candidates
            current_frontier = top_w_count_kmers<std::vector<unsigned int>, float>(new_frontier, width);

            if (DEBUG) {
                seqan3::debug_stream << "Frontier: ";
                for (auto [path, weight] : current_frontier) {
                    seqan3::debug_stream << kmer_vector_to_sequence(path, k) << " " << weight << "\n";
                }
                seqan3::debug_stream << "\n";
            }
        }

        // the best path is simply the longest path with the highest weight
        float max_weight = 0;
        unsigned int max_length = 0;
        bool found_target = false;
        std::vector<unsigned int> max_path;
        for (auto [path, weight] : current_frontier) {
            //seqan3::debug_stream << "Path: " << kmer_vector_to_sequence(path, k) << " " << weight << "\n";
            //seqan3::debug_stream << kmer_vector_to_sequence(max_path, k) << ", " <<  max_weight << "\n";
            bool match_target = (path[path.size() - 1] == target_kmer);
            if ((!found_target && match_target) || (((found_target && match_target) || (!found_target && !match_target)) && weight > max_weight)) {
                max_weight = weight;
                max_length = path.size();
                max_path = path;
                if (match_target) {
                    found_target = true;
                }
            }
        }
        //seqan3::debug_stream << "Max kmer" << hash_to_kmer(max_kmer, k) << " " << max_weight << "\n";

        // Turn the path into a sequence
        seqan3::dna4_vector res = kmer_vector_to_sequence(max_path, k);
        //seqan3::debug_stream << "RES:" << res << "\n";


        return res;
    }



    seqan3::dna4_vector find_consensus_beam_bidirectional(const std::vector<seqan3::dna4_vector>& sv, uint8_t k, unsigned int width, unsigned int expected_length) {
        // find the initial k-mer 
        auto [initial_counts, last_counts] = initial_and_last_kmer_counts(sv, k);
        auto all_counts = all_kmer_counts(sv, k);

        auto number_of_sequences = sv.size();

        // greedily find the next k-mer.
        auto initial_kmer = max_count_kmer(initial_counts);
        auto target_kmer = max_count_kmer(last_counts);

        auto top_w_init_kmers = top_w_count_kmers<unsigned int, unsigned int>(initial_counts, width);

        // initialize variables that store the path to the k-mer
        std::vector<std::pair<std::vector<unsigned int>, float>> forward_frontier;

        for (auto const [kmer, weight] : top_w_init_kmers) {
            forward_frontier.push_back(std::make_pair(std::vector<unsigned int>{kmer}, std::log2((float)weight)));
        }

        //forward_frontier.push_back(std::make_pair(std::vector<unsigned int>{initial_kmer}, std::log2((float)initial_counts[initial_kmer])));

        // forward beam search
        unsigned int forward_length = (expected_length - k + 1) / 2;
        unsigned int backward_length = expected_length - k + 1 - forward_length;
        //seqan3::debug_stream << "Forward length: " << forward_length << "\n";
        //seqan3::debug_stream << "Backward length: " << backward_length << "\n";



        // Perform greedy beam search
        for (unsigned int l = 0; l < forward_length; l++) {
            auto new_frontier = beam_search_step(forward_frontier, all_counts, width, true, k);

            if (new_frontier.empty()) {
                // no more candidates
                break;
            }
            forward_frontier = new_frontier;
        }


        // backward beam search
        std::vector<std::pair<std::vector<unsigned int>, float>> backward_frontier;

        auto top_w_target_kmers = top_w_count_kmers<unsigned int, unsigned int>(last_counts, width);


        for (auto const [kmer, weight] : top_w_target_kmers) {
            backward_frontier.push_back(std::make_pair(std::vector<unsigned int>{kmer}, std::log2((float)weight)));
        }

        //backward_frontier.push_back(std::make_pair(std::vector<unsigned int>{target_kmer}, std::log2((float)last_counts[target_kmer])));

        for (unsigned int l = 0; l < backward_length - 1; l++) {
            auto new_frontier = beam_search_step(backward_frontier, all_counts, width, false, k);

            if (new_frontier.empty()) {
                // no more candidates
                break;
            }
            backward_frontier = new_frontier;
        }

        // merge the two frontiers
        // create a map that contains the common k-mers between the two frontiers
        std::unordered_map<unsigned int, float> forward_frontier_kmer_to_weight;
        std::unordered_map<unsigned int, std::vector<unsigned int>> forward_frontier_kmer_to_path;

        for (auto [path, weight] : forward_frontier) {
            auto kmer = path[path.size() - 1];
            if (forward_frontier_kmer_to_weight.contains(kmer)) {
                if (weight > forward_frontier_kmer_to_weight[kmer]) {
                    forward_frontier_kmer_to_weight[kmer] = weight;
                    forward_frontier_kmer_to_path[kmer] = path;
                }
            } else {
                forward_frontier_kmer_to_weight[kmer] = weight;
                forward_frontier_kmer_to_path[kmer] = path;
            }
        }

        // check the backward frontier
        std::unordered_map<unsigned int, float> backward_frontier_kmer_to_weight;
        std::unordered_map<unsigned int, std::vector<unsigned int>> backward_frontier_kmer_to_path;

        for (auto [path, weight] : backward_frontier) {
            auto kmer = path[path.size() - 1];
            if (backward_frontier_kmer_to_weight.contains(kmer)) {
                if (weight > backward_frontier_kmer_to_weight[kmer]) {
                    backward_frontier_kmer_to_weight[kmer] = weight;
                    backward_frontier_kmer_to_path[kmer] = path;
                }
            } else {
                backward_frontier_kmer_to_weight[kmer] = weight;
                backward_frontier_kmer_to_path[kmer] = path;
            }
        }

        // merge the two frontiers
        // pick the k-mer with the highest weight
        float max_weight = 0;
        unsigned int max_length = 0;
        bool found_target = false;
        std::vector<unsigned int> max_path;
        for (auto [kmer, weight] : forward_frontier_kmer_to_weight) {
            if (backward_frontier_kmer_to_weight.contains(kmer)) {
                float total_weight = weight + backward_frontier_kmer_to_weight[kmer];
                if (total_weight > max_weight) {
                    max_weight = total_weight;
                    max_path = forward_frontier_kmer_to_path[kmer];
                    auto& backward_path = backward_frontier_kmer_to_path[kmer];
                    max_path.insert(max_path.end(), backward_path.rbegin() + 1, backward_path.rend());
                }
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
        //return find_consensus_beam(sv, k, 50, 108);
        return find_consensus_beam_bidirectional(sv, k, 50, 108);
        //return find_consensus_bfs(sv, k, 110);
        //return find_consensus_greedy(sv, k, 150);
    }
};

#endif