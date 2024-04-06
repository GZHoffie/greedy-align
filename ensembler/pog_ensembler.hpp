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
    bool is_start;               // whether this node is a start node

    // list of neighbors, with the second element being the edge weight 
    std::unordered_map<std::shared_ptr<po_node>, unsigned int> next;  

    po_node(unsigned int kmer_length, unsigned int position, unsigned int kmer_key, unsigned int consecutive_occurance = 1) {
        k = kmer_length;
        offset = position;
        key = kmer_key;
        repetition = consecutive_occurance;
        is_start = false;
    }

    po_node() : po_node(0, 0, 0, 0) {
        is_start = true;
    }

    ~po_node() = default;

    void print() {
        seqan3::debug_stream << "Offset: " << offset << ",\tK-mer: " << hash_to_kmer(key, k) << ",\tNumber of repetitions: " << repetition << ",\tNeighbors: ";
        for (auto& nl : next) {
            seqan3::debug_stream << hash_to_kmer(nl.first->key, k) << " (" << nl.second << "); ";
        }
        seqan3::debug_stream << "\n";
    }
};

class po_graph {
private:
    unsigned int k;
    unsigned int error_tolerence;
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
    std::shared_ptr<po_node> _add_new_node(std::shared_ptr<po_node> last_node, unsigned int position, unsigned int kmer_key, unsigned int consecutive_occurence = 1) {
        // STEP 1: check if this node is a neighbor of the last node
        // if yes, return, and no new nodes are created
        if (last_node != nullptr) {
            for (auto & [neighbor, weight] : last_node->next) {
                if (neighbor->key == kmer_key && neighbor->repetition == consecutive_occurence) {
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
                if (n->repetition == consecutive_occurence && std::abs((int) n->offset - (int) position) <= min_distance) {
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
        auto new_node = std::make_shared<po_node>(k, position, kmer_key, consecutive_occurence);
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
        for (unsigned int i = 0; i < highways.size(); i++) {
            // the position of this highway
            unsigned int highway_offset = highways[i].offset;

            // get the k-mers in the highway
            auto kmer_hash = consensus[i] | seqan3::views::kmer_hash(seqan3::ungapped{k});
            std::vector<unsigned int> kmer_hash_vector(kmer_hash.begin(), kmer_hash.end());

            unsigned int j = 0;
            std::shared_ptr<po_node> last_node = nullptr;
            while (j < kmer_hash_vector.size()) {
                // offset of this k-mer
                unsigned int kmer_offset = highway_offset + j;

                // key of this k-mer
                unsigned int kmer_key = kmer_hash_vector[j];

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
                last_node = _add_new_node(last_node, kmer_offset, kmer_key, repetition);
                j++;
            }
        }
    }


    seqan3::dna4_vector _path_to_dna4_vector(const std::vector<unsigned int>& kmers) {
        using namespace seqan3::literals;

        seqan3::dna4_vector res;
        if (kmers.empty()) return res;
        
        res = hash_to_kmer(kmers[0], k);
        for (unsigned int i = 1; i < kmers.size(); i++) {
            auto rank = kmers[i] & 3;
                switch (rank){
                case 0:
                    res.push_back('A'_dna4);
                    break;
                case 1:
                    res.push_back('C'_dna4);
                    break;
                case 2:
                    res.push_back('G'_dna4);
                    break;
                default:
                    res.push_back('T'_dna4);
                    break;
            }
        }
        return res;
    }

public:
    po_graph(unsigned int kmer_length, unsigned int max_error) {
        k = kmer_length;
        error_tolerence = max_error;
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
     */
    seqan3::dna4_vector consensus_path() {
        std::vector<unsigned int> res;
        auto current_node = start;
        while (!current_node->next.empty()) {
            unsigned int highest_weight = 0;
            std::shared_ptr<po_node> best_neighbor;
            for (auto neighbor : current_node->next) {
                if (neighbor.second > highest_weight) {
                    best_neighbor = neighbor.first;
                }
            }
            for (unsigned int i = 0; i < best_neighbor->repetition; i++) {
                res.push_back(best_neighbor->key);
            }
            current_node = best_neighbor;
        }
        return _path_to_dna4_vector(res);
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

    consensus_t _highway_to_concensus(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2, const std::vector<de_bruijn_highway_t> highways) {
        // vectors to store the concensus
        std::vector<seqan3::dna4_vector> c1, c2;
        
        for (auto & h : highways) {
            seqan3::dna4_vector h1(&s1[h.offset - h.lane], &s1[h.offset - h.lane + h.length + align_k - 1]);
            seqan3::dna4_vector h2(&s2[h.offset], &s2[h.offset + h.length + align_k - 1]);

            seqan3::debug_stream << h1 << "\n" << h2 << "\n";
            c1.push_back(h1);
            c2.push_back(h2);
        }

        return std::make_pair(c1, c2);
    }


public:
    partial_order_graph_ensembler(unsigned int pog_kmer_length, unsigned int align_kmer_length, unsigned int band_width, 
                                  unsigned int max_errors, bool print_debug_messages = false) {
        k = pog_kmer_length;
        align_k = align_kmer_length;
        aligner = new greedy_aligner<READ_LENGTH>(align_kmer_length, band_width, max_errors, print_debug_messages);
        graph = new po_graph(pog_kmer_length, band_width);
    }

    ~partial_order_graph_ensembler() {
        delete aligner;
        delete graph;
    }

    seqan3::dna4_vector ensemble(const std::vector<seqan3::dna4_vector>& sv) {
        // TODO: try sampling pairs of sequences for ensembling.
        seqan3::debug_stream << sv.size() << "\n";
        int num_sequences = sv.size();

        unsigned i = 0, j = 1;
        auto s1 = sv[i];
        auto s2 = sv[j];

        // find consensus between two highways
        auto highways = aligner->find_long_consecutive_matches(s1, s2);
        auto consensus = _highway_to_concensus(s1, s2, highways);
        seqan3::debug_stream << "OK\n";

        graph->read_consensus(highways, consensus);
        graph->print();

        auto res = graph->consensus_path();
        seqan3::debug_stream << res << "\n";


        //seqan3::dna4_vector res;
        return res;
    }

    


    

};

#endif