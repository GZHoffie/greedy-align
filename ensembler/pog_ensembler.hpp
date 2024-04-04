#ifndef GREEDY_ALIGN_POG_ENSEMBLER_HPP
#define GREEDY_ALIGN_POG_ENSEMBLER_HPP


#include "ensembler.hpp"
#include "../aligner/greedy_aligner.hpp"

#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <bitset>

typedef std::pair<std::vector<seqan3::dna4_vector>, std::vector<seqan3::dna4_vector>> consensus_t;

class po_node {
private:
    unsigned int offset;
    unsigned int key;
    unsigned int index;
    bool is_start;
    std::vector<po_node*> next;

public:
    po_node(unsigned int position, unsigned int kmer_key, unsigned int occurance = 0) {
        offset = position;
        key = kmer_key;
        index = occurance;
        is_start = false;
    }

    po_node() : po_node(0, 0, 0) {
        is_start = true;
    }
};

class po_graph {
private:
    unsigned int k;
    po_node start;

    // map that allow quick search of nodes.
    std::unordered_map<unsigned int, std::vector<po_node*>> nodes;

public:
    void read_consensus(const consensus_t& consensus) {
        auto &[c1, c2] = consensus;
        std::unordered_map<unsigned int, unsigned int> kmer_occurences;
        po_node* prev = &start; 

        for (auto &c : c1) {
            auto kmer_keys = c | seqan3::views::kmer_hash(seqan3::ungapped{k});
            for (auto kmer_key: kmer_keys) {
                // check how many time this k-mer has occurred.
                unsigned int occurence = 0;
                if (kmer_occurences.contains(kmer_key)) {
                    occurence = kmer_occurences[kmer_key];
                } 
                kmer_occurences[kmer_key] = occurence + 1;

                // check if it matches one of the candidates

            }

        }
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
    }

    ~partial_order_graph_ensembler() {
        delete aligner;
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

        seqan3::dna4_vector res;
        return res;
    }


    

};

#endif