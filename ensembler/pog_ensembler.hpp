#ifndef GREEDY_ALIGN_POG_ENSEMBLER_HPP
#define GREEDY_ALIGN_POG_ENSEMBLER_HPP


#include "ensembler.hpp"
#include "../aligner/greedy_aligner.hpp"

#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <bitset>

typedef std::pair<std::vector<seqan3::dna4_vector>, std::vector<seqan3::dna4_vector>> consensus_t;

class po_node {
public:
    // node attributes
    unsigned int offset;         // the rough position of the k-mer in the traces
    unsigned int key;            // an integer (hash value) representing which k-mer the node is representing, using seqan3::kmer_hash
    unsigned int repetition;     // how many time this k-mer occurs consecutively. For homopolymers
    unsigned int k;              // length of the k-mer
    unsigned int occurence;      // number of time this k-mer has occurred in the string
    bool is_start;               // whether this node is a start node
    bool visited;                // whether the node is visited

    // list of neighbors, with the second element being the edge weight 
    std::unordered_map<std::shared_ptr<po_node>, unsigned int> next;  

    po_node(unsigned int kmer_length, unsigned int position, unsigned int kmer_key, unsigned int consecutive_occurence = 1, unsigned int total_occurence = 0) {
        k = kmer_length;
        offset = position;
        key = kmer_key;
        repetition = consecutive_occurence;
        occurence = total_occurence;
        is_start = false;
        visited = false;
    }

    po_node() : po_node(0, 0, 0, 0, 0) {
        is_start = true;
        visited = false;
    }

    ~po_node() {
        //seqan3::debug_stream << "Destroyed\n";
    }

    void print() {
        seqan3::debug_stream << "Offset: " << offset << ",\tK-mer: " << hash_to_kmer(key, k) << ",\trepetitions: " << repetition << ",\toccurences: " << occurence << ",\t visited:" << visited << ",\tNeighbors: ";
        for (auto& nl : next) {
            seqan3::debug_stream << hash_to_kmer(nl.first->key, k) << "-" << nl.first->offset << "," << nl.first->occurence  << "," << nl.first->repetition << " (" << nl.second << "); ";
        }
        seqan3::debug_stream << "\n";
    }
};

class po_graph {
private:
    unsigned int k;
    unsigned int error_tolerence;
    unsigned int target_length;
    std::shared_ptr<po_node> start;

    // map that allow quick search of nodes.
    std::unordered_map<unsigned int, std::vector<std::shared_ptr<po_node>>> nodes_list;

    /**
     * @brief Add a new node to the graph.
     * If the node already exists in the graph, then regard that node as the to-be-added new node.
     * Increase the edge weight if last_node is specified.
     * 
     * @param last_node the pointer to the last node, if it exists. Otherwise, it should be a nullptr.
     * @param position the rough position of this node in the consensus.
     * @param kmer_key the hash value of the kmer the node is representing.
     * @param consecutive_occurence how many times the k-mer repeats (for homopolymers)
     */
    std::shared_ptr<po_node> _add_new_node(std::shared_ptr<po_node> last_node, unsigned int position, unsigned int kmer_key, unsigned int consecutive_occurence = 1, unsigned int total_occurence = 0) {
        // STEP 1: check if this node is a neighbor of the last node
        // if yes, return, and no new nodes are created
        if (last_node != nullptr) {
            for (auto & [neighbor, weight] : last_node->next) {
                if (neighbor->key == kmer_key && 
                    neighbor->repetition == consecutive_occurence &&
                    neighbor->occurence == total_occurence) {
                    last_node->next[neighbor] += 1;
                    return neighbor;
                }
            }
        }

        // STEP 2: check if the new node is similar to one of the previously created nodes.
        if (nodes_list.contains(kmer_key)) {
            std::shared_ptr<po_node> closest_node = nullptr;
            unsigned int min_distance = error_tolerence;
            for (auto &n : nodes_list[kmer_key]) {
                if (n->repetition == consecutive_occurence && 
                    n->occurence == total_occurence &&
                    std::abs((int) n->offset - (int) position) <= min_distance) {
                    closest_node = n;
                }
            }
            // if found a very similar node, simply regard the node we are adding is this node.
            if (closest_node != nullptr) {
                // add this to the neighbor of last_node
                if (last_node != nullptr) {
                    last_node->next[closest_node] = 1;
                }
                return closest_node;
            }
        }

        // STEP 3: we add a new node if the above cases don't happen.
        auto new_node = std::make_shared<po_node>(k, position, kmer_key, consecutive_occurence, total_occurence);
        // add this to the list of neighbors of last_node
        if (last_node != nullptr) {
            last_node->next[new_node] = 1;
        }
        // add the new node to the node list
        if (nodes_list.contains(kmer_key)) {
            nodes_list[kmer_key].push_back(new_node);
        } else {
            std::vector<std::shared_ptr<po_node>> l{new_node};
            nodes_list[kmer_key] = l;
        }
        return new_node;
    }

    void _add_consensus(const std::vector<de_bruijn_highway_t>& highways, const std::vector<seqan3::dna4_vector>& consensus) {
        assert(highways.size() == consensus.size());

        // Map that stores the occurence of k-mers
        std::unordered_map<unsigned int, unsigned int> occurences;


        for (unsigned int i = 0; i < highways.size(); i++) {
            // the position of this highway
            unsigned int highway_offset = highways[i].offset;

            // get the k-mers in the highway
            auto kmer_hash = consensus[i] | seqan3::views::kmer_hash(seqan3::ungapped{(uint8_t)k});
            std::vector<unsigned int> kmer_hash_vector(kmer_hash.begin(), kmer_hash.end());

            unsigned int j = 0;
            std::shared_ptr<po_node> last_node = nullptr;
            while (j < kmer_hash_vector.size()) {
                // offset of this k-mer
                unsigned int kmer_offset = highway_offset + j;

                // key of this k-mer
                unsigned int kmer_key = kmer_hash_vector[j];

                // number of time this k-mer has occurred
                unsigned int kmer_occurence = 0;
                if (occurences.contains(kmer_key)) {
                    kmer_occurence = occurences[kmer_key] + 1;
                }
                occurences[kmer_key] = kmer_occurence;

                // check if this is a homopolymer
                unsigned int repetition = 1;
                while (j + 1 < kmer_hash_vector.size() && kmer_hash_vector[j + 1] == kmer_key) {
                    repetition++;
                    j++;
                }

                // insert the new node.
                if (kmer_offset == 0) {
                    last_node = start;
                }
                last_node = _add_new_node(last_node, kmer_offset, kmer_key, repetition, kmer_occurence);
                j++;
            }
        }
    }

    /**
     * Clean the partial order graph.
     * * Merge similar nodes that has the common ancestor together.
     * * Clean the edges that have small weight.
     * 
     * @param min_votes minimum number of votes to regard that edge as valid
     */
    void _clean_graph(unsigned int min_votes) {
        // TODO: STEP 1: Merge nodes together

        // STEP 2: Remove the edges that have less weight than min_votes
        for (auto & item : nodes_list) {
            for (auto & node : item.second) {
                // remove neighbor if it doesn't have enough weight
                auto begin_iter = node->next.begin();
                auto end_iter = node->next.end();
                for (auto iter = begin_iter; iter != end_iter;) {
                    if (iter->second < min_votes) {
                        iter = node->next.erase(iter);
                    } else {
                        ++iter;
                    }
                }
            }
        }
        // for the start node also
        auto begin_iter = start->next.begin();
        auto end_iter = start->next.end();
        for (auto iter = begin_iter; iter != end_iter;) {
            if (iter->second < min_votes) {
                iter = start->next.erase(iter);
            } else {
                ++iter;
            }
        }
    }

    typedef std::pair<std::vector<std::shared_ptr<po_node>>, unsigned int> path_t;

    std::vector<path_t> all_paths;

    void _dfs_util(std::shared_ptr<po_node> current_node, std::vector<std::shared_ptr<po_node>>& current_path, unsigned int current_path_votes) {
        //seqan3::debug_stream << current_path << " " << current_path_votes << "\n";
        // start the recursive step
        if (current_node->next.empty()) {
            // end of the path, store the path in all_paths
            all_paths.push_back(std::make_pair(current_path, current_path_votes));
        } else {
            // if current_path already contains the current node, this path is invalid.
            //seqan3::debug_stream << current_path << " "  << "\n";
            //current_node->print();
            if (std::find(current_path.begin(), current_path.end() - 1, current_node) != current_path.end() - 1) return;

            // otherwise, continue the search
            for (auto& item : current_node->next) {
                // Add this node to the path
                std::vector<std::shared_ptr<po_node>> path(current_path);
                unsigned int path_votes = current_path_votes + item.second;
                path.push_back(item.first);
                //seqan3::debug_stream << "After: " << path << " " << path_votes << " " <<  item.first << "\n";

                // Do DFS after that node
                _dfs_util(item.first, path, path_votes);
            }
        }
    }

    void _dfs() {
        // reset contents of the paths
        all_paths.clear();


        // start the recursion
        std::vector<std::shared_ptr<po_node>> path;
        unsigned int path_votes = 0;
        _dfs_util(start, path, path_votes);

        //seqan3::debug_stream << all_paths << "\n";
    }


    seqan3::dna4_vector _path_to_dna4_vector(const std::vector<std::shared_ptr<po_node>>& path) {
        using namespace seqan3::literals;

        seqan3::dna4_vector res;
        if (path.empty()) return res;
        unsigned int last_hash; // hash function of the last seen k-mer
        
        last_hash = path[0]->key;
        res = hash_to_kmer(path[0]->key, k);
        for (unsigned int i = 1; i < path.size(); i++) {
            auto rank = path[i]->key & 3;
            auto nt = seqan3::assign_rank_to(rank, seqan3::dna4{});
            for (unsigned int j = 0; j < path[i]->repetition; j++)
                res.push_back(nt);
        }
        //seqan3::debug_stream << res << "\n";
        return res;
    }


    seqan3::dna4_vector _find_best_path() {
        // storing the best path
        int min_length_diff = target_length;
        unsigned int max_votes = 0;
        seqan3::dna4_vector res;

        // iterate through all possible solutions
        for (auto &[path, path_weight] : all_paths) {
            auto dna4_path = _path_to_dna4_vector(path);
            int length_diff = std::abs((int) dna4_path.size() - (int)target_length);
            if (length_diff < min_length_diff || (length_diff == min_length_diff && path_weight > max_votes)) {
                res = dna4_path;
                min_length_diff = length_diff;
                max_votes = path_weight;
            }
        }

        return res;
    }
public:
    po_graph(unsigned int kmer_length, unsigned int max_error, unsigned int consensus_length) {
        k = kmer_length;
        error_tolerence = max_error;
        target_length = consensus_length;
        start = std::make_shared<po_node>();
    }

    void read_consensus(const std::vector<de_bruijn_highway_t>& highways, const consensus_t& consensus) {
        auto &[c1, c2] = consensus;
        _add_consensus(highways, c1);
        _add_consensus(highways, c2);
    }

    void print() {
        seqan3::debug_stream << "Nodes:\n";
        for (auto nl : nodes_list) {
            seqan3::debug_stream << "- K-Mer: " << hash_to_kmer(nl.first, k) << "\n";
            for (auto n : nl.second) {
                n->print();
            }
        }
    }

    /**
     * @brief find the path that has the highest weight.
     * Use DFS to find a valid path, that is closest to the target length
     */
    seqan3::dna4_vector consensus_path(unsigned int num_sequences) {
        std::vector<unsigned int> res;

        // TODO: find a threshold to recognize a correct path. Here this number is just random:
        // At least 10% of all the pairs (n*(n-1)/2) agree with the path.
        unsigned int min_path_weight = (unsigned int) (num_sequences * (num_sequences - 1) / 20);
        //this->print();
        _clean_graph(min_path_weight);

        //this->print();

        // perform DFS
        _dfs();

        return _find_best_path();
    }
    
    void reset() {
        start = std::make_shared<po_node>();
        // clean all references
        // TODO: check if all references are cleaned.
        for (auto & item: nodes_list) {
            for (auto & node: item.second) {
                node->next.clear();
            }
        }
        nodes_list.clear();
        //seqan3::debug_stream << "==========================\n";
    }


};

template<unsigned int READ_LENGTH>
class partial_order_graph_ensembler : public ensembler {
private:
    // k-mer length used to construct the partial order graph
    unsigned int k; 
    unsigned int align_k;

    // pairwise alignment to find long consecutive matches
    greedy_aligner<READ_LENGTH>* aligner;

    // partial order graph to ensemble the consensus
    po_graph* graph;

    // target length of the consensus
    unsigned int target_length;

    bool debug;

    consensus_t _highway_to_concensus(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2, const std::vector<de_bruijn_highway_t> highways) {
        // vectors to store the concensus
        std::vector<seqan3::dna4_vector> c1, c2;
        
        for (auto & h : highways) {
            seqan3::dna4_vector h1(&s1[std::max(0, (int)h.offset - (int)h.lane)], &s1[h.offset - h.lane + h.length + align_k - 1]);
            seqan3::dna4_vector h2(&s2[h.offset], &s2[h.offset + h.length + align_k - 1]);

            if (debug) {
                seqan3::debug_stream << h1 << "\n" << h2 << "\n";
            }
            c1.push_back(h1);
            c2.push_back(h2);
        }

        return std::make_pair(c1, c2);
    }


public:
    partial_order_graph_ensembler(unsigned int pog_kmer_length, unsigned int align_kmer_length, unsigned int band_width, 
                                  unsigned int max_errors, unsigned int target_consensus_length, bool print_debug_messages = false) {
        k = pog_kmer_length;
        align_k = align_kmer_length;
        aligner = new greedy_aligner<READ_LENGTH>(align_kmer_length, band_width, max_errors, print_debug_messages);
        graph = new po_graph(pog_kmer_length, band_width, target_consensus_length);
        debug = print_debug_messages;
    }

    ~partial_order_graph_ensembler() {
        delete aligner;
        delete graph;
    }

    seqan3::dna4_vector ensemble(const std::vector<seqan3::dna4_vector>& sv) {
        graph->reset();
        // TODO: try sampling pairs of sequences for ensembling.
        for (unsigned int i = 0; i < sv.size(); i++) {
            for (unsigned int j = i + 1; j < sv.size(); j++) {
                //seqan3::debug_stream << i << " " << j << "\n";
                auto s1 = sv[i];
                auto s2 = sv[j];

                // find consensus between two highways
                auto highways = aligner->find_long_consecutive_matches(s1, s2);
                auto consensus = _highway_to_concensus(s1, s2, highways);
        
                graph->read_consensus(highways, consensus);
            }
        }
        //graph->print();
        auto res = graph->consensus_path(sv.size());
        //seqan3::debug_stream << res << "\n";
        if (debug) {
            graph->print();
            seqan3::debug_stream << res << "\n";
        }


        //seqan3::dna4_vector res;
        return res;
    }

    


    

};

#endif