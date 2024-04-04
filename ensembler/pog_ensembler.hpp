#ifndef GREEDY_ALIGN_POG_ENSEMBLER_HPP
#define GREEDY_ALIGN_POG_ENSEMBLER_HPP


#include "ensembler.hpp"
#include "../aligner/greedy_aligner.hpp"

#include <seqan3/search/kmer_index/shape.hpp>
#include <bitset>

template<unsigned int READ_LENGTH>
class partial_order_graph_ensembler : public ensembler {
private:
    // k-mer length used to construct the partial order graph
    unsigned int k; 
    unsigned int align_k;

    // pairwise alignment to find long consecutive matches
    greedy_aligner<READ_LENGTH>* aligner;

    std::pair<std::vector<seqan3::dna4_vector>, std::vector<seqan3::dna4_vector>> 
    _highway_to_concensus(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2, const std::vector<de_bruijn_highway_t> highways) {
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