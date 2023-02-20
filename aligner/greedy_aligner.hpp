#ifndef GREEDY_ALIGN_GREEDY_ALIGNER_HPP
#define GREEDY_ALIGN_GREEDY_ALIGNER_HPP


#include "aligner.hpp"

#include <seqan3/search/kmer_index/shape.hpp>
#include <bitset>

using seqan3::operator""_cigar_operation;


template<unsigned int READ_LENGTH>
class read_bit_vectors {
public:
    std::bitset<READ_LENGTH> b1;
    std::bitset<READ_LENGTH> b2;
    unsigned int length;

    /**
     * @brief Construct a new read bit vectors object using a DNA sequence.
     * 
     * @param s the seqan3::dna4_vector object, representing the sequence to
     *          be converted to bit representation.
     */
    read_bit_vectors(const seqan3::dna4_vector& s) {
        for (int i = 0; i < s.size(); i++) {
            b1[i] = s[i].to_rank() & 1;
            b2[i] = s[i].to_rank() & 2;
        }
        length = s.size();
        
    }

    std::bitset<READ_LENGTH> shifted_xor(const read_bit_vectors<READ_LENGTH>& that, int shift_amount) {
        std::bitset<READ_LENGTH> res;
        // void out all meaningless bits
        res = ~res << (READ_LENGTH - that.length);
        //seqan3::debug_stream << shift_amount << " " << res << "\n";

        if (shift_amount > 0) {
            // void out the shifted bits
            res = ~(res >> (READ_LENGTH - that.length) >> shift_amount << shift_amount);
            //seqan3::debug_stream << shift_amount << " " << res << "\n";
            // perform bitwise xor operation
            res |= ((b1 << shift_amount) ^ that.b1) | ((b2 << shift_amount) ^ that.b2);
            //seqan3::debug_stream << shift_amount << " " << res << "\n";
        } else {
            // same as above, but in opposite direction
            res = ~(res << -shift_amount >> (-shift_amount + READ_LENGTH - that.length));
            //seqan3::debug_stream << shift_amount << " " << res << "\n";
            res |= ((b1 >> -shift_amount) ^ that.b1) | ((b2 >> -shift_amount) ^ that.b2);
            //seqan3::debug_stream << shift_amount << " " << res << "\n";
        }
        return res;
    }

};


typedef struct {
    int lane;
    unsigned int offset;
    unsigned int length = 0;
} de_bruijn_highway_t;


template<unsigned int READ_LENGTH>
class de_bruijn_lanes {
private:
    // bit vectors storing the highways
    std::vector<std::bitset<READ_LENGTH>> lanes;
    std::vector<std::bitset<READ_LENGTH>> lanes_flipped;

    // bandwidth, the maximum number of indels allowed.
    int bw;

    // k, the length of k-mers used to match the two sequences.
    unsigned int k;


    /**
     * @brief After assigning all the bits in `lanes`, store the ~`lanes`
     *        in `lanes_flipped`, so that searching for highways is faster.
     */
    void _store_flipped_lanes() {
        for (int i = 0; i < lanes.size(); i++) {
            lanes_flipped[i] = ~lanes[i];
        }
    }

    /**
     * @brief Check the hamming distance between (v1 << shift_amount) ^ v2, and store
     *        the results in `lanes`.
     * 
     * @param v1, v2 the two sequences, already converted to bit vectors.
     * @param shift_amount The amount of shifting, which should be in [-band_width, band_width].
     */
    void _shifted_xor(read_bit_vectors<READ_LENGTH>& v1, 
                      read_bit_vectors<READ_LENGTH>& v2, 
                      int shift_amount) {
    
        lanes[shift_amount + bw] = v1.shifted_xor(v2, shift_amount);
    }


public:
    // length of the read and pattern,
    unsigned int l1, l2;

    /**
     * @brief Construct a new de_bruijn_lanes object.
     * 
     * @param band_width Maximum number of indels allowed.
     * @param seed_length length of k-mers used in de Bruijn path.
     */
    de_bruijn_lanes(unsigned int seed_length, int band_width) {
        bw = band_width;
        for (int i = -band_width; i <= band_width; i++) {
            std::bitset<READ_LENGTH> lane, lane_rev;
            lanes.push_back(lane);
            lanes_flipped.push_back(lane_rev);
        }
        k = seed_length;
    }

    /**
     * @brief Clear all content in the bit vectors.
     */
    void clear() {
        for (unsigned int i = 0; i < lanes.size(); i++) {
            lanes[i].reset();
        }
    }

    /**
     * @brief Get the corresponding bitset with shift amount i.
     * 
     * @param i a value between [-band_width, band_width]
     * @return std::bitset<READ_LENGTH>& the reference to the bitset.
     */
    std::bitset<READ_LENGTH>& get(int i) {
        return lanes[i + bw];
    }

    /**
     * @brief Print out the lanes in order, for debugging.
     */
    void print() {
        seqan3::debug_stream << "[INFO]\t\t" << "Current de Bruijn lane layout:\n";
        for (int i = -bw; i <= bw; i++) {
            seqan3::debug_stream << "[INFO]\t\t" << "Lane " << i << "\t| " << get(i) << "\n";
        }
    }


    /**
     * @brief Find highways (common subpaths) in the de Bruijn graphs, which are represented by
     *        (consecutive) zeros in the lanes.
     * 
     * @param s1, s2 the DNA sequences represented using seqan3::dna4 objects.
     */
    void read(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) {
        // build read_bit_vectors objects base on the reads
        read_bit_vectors<READ_LENGTH> v1(s1);
        read_bit_vectors<READ_LENGTH> v2(s2);

        // store the length of the two sequences.
        l1 = s1.size() - k + 1;
        l2 = s2.size() - k + 1;

        // perform shifting and bitwise xor operation, and store in `lanes`.
        for (int l = -bw; l <= bw; l++) {
            _shifted_xor(v1, v2, l);
            auto original_lane = get(l);
            for (int i = 1; i < k; i++) {
                get(l) |= (original_lane >> i);
            }
        }

        // store flipped bits in `lanes_flipped`.
        _store_flipped_lanes();
    }


    
    /**
     * @brief Find the position and length of the next highway in lane after offset.
     * 
     * @param lane the index for the lane.
     * @param offset the position where we start to look for a highway.
     * @return the information about the starting position and the second represent the length of highway.
     */
    de_bruijn_highway_t nearest_highway(int lane, unsigned int offset) {
        de_bruijn_highway_t res;
        res.lane = lane;

        // find the starting point of the highway
        if (offset == 0) {
            res.offset = lanes_flipped[lane + bw]._Find_first();
        } else {
            res.offset = lanes_flipped[lane + bw]._Find_next(offset - 1);
        }
        
        // find the length of the highway
        if (res.offset < READ_LENGTH) res.length = lanes[lane + bw]._Find_next(res.offset) - res.offset;
        else res.length = 0;

        seqan3::debug_stream << "[INFO]\t\t" << "On lane " << res.lane << ", next highway starts at " 
                             << res.offset << " and has length " << res.length << ".\n";
        return res;
    }
};

template<unsigned int READ_LENGTH>
class greedy_aligner : public aligner {
private:
    unsigned int k;
    int bw;
    unsigned int e;
    de_bruijn_lanes<READ_LENGTH>* lanes;

    /**
     * @brief calculate the offset if we go from current lane to a target lane.
     * 
     * @param current_lane the current de Bruijn lane we are in.
     * @param target_lane the target lane we are going to (by insertion or deletion).
     * @param current_offset the current k-mer we are at.
     */
    unsigned int _get_offset(int current_lane, int target_lane, unsigned int current_offset) {
        // if at the beginning of the alignment, set offset in all lanes to zero.
        if (current_offset == 0) return 0;
        // if there is no insertion/deletion, make offset the same as before.
        if (current_lane == target_lane) return current_offset;
        // if there are insertion/deletions, increase the offset by k-1+|current_lane - target_lane|.
        else return current_offset + k + std::abs(current_lane - target_lane) - 1;
    }

    void _update_cigar(CIGAR& cigar, int current_lane, unsigned int current_offset, const de_bruijn_highway_t& best_highway) {
        // record number of indels
        if (current_lane > best_highway.lane) {
            cigar.push_back(current_lane - best_highway.lane, 'D'_cigar_operation);
        } else if (current_lane < best_highway.lane) {
            cigar.push_back(best_highway.lane - current_lane, 'I'_cigar_operation);
        }
        
        // find number of mismatches
        int mismatches = best_highway.offset - current_offset - std::abs(current_lane - best_highway.lane) - k + 1;
        if (mismatches > 0) {
            cigar.push_back(mismatches, 'X'_cigar_operation);
        }
        cigar.push_back(best_highway.length, '='_cigar_operation);
    }


    /**
     * @brief Greedily find the subpath that maximize the k-matching.
     * 
     */
    align_result_t _find_subpath() {
        // initialize result
        align_result_t res;
        CIGAR cigar;

        // record the current position
        int current_lane = 0, current_offset = 0;

        // greedily search for the next best highway
        unsigned int best_errors;
        while (true) {
            // initialize
            best_errors = INT_MAX;
            de_bruijn_highway_t best_highway;

            // check each lane for possible lengths
            for (int l = -bw; l <= bw; l++) {
                auto offset_l = _get_offset(current_lane, l, current_offset);
                auto highway_l = lanes->nearest_highway(l, offset_l);

                // skip if there is no highway in this lane
                if (highway_l.length == 0) continue;

                // the number of k-mers we missed
                unsigned int errors_l = highway_l.offset - current_offset;

                // chech if the current highway is a better one.
                if (errors_l <= std::max(e, best_errors) && 
                    (best_errors > e || 
                     highway_l.length > best_highway.length ||
                     (highway_l.length == best_highway.length && errors_l < best_errors)
                    )) {
                    best_highway = highway_l;
                    best_errors = errors_l;
                }
                
            }
            seqan3::debug_stream << "[INFO]\t\tChosen next highway to be on lane " << best_highway.lane << " starting at " << best_highway.offset 
                                 << " with length " << best_highway.length << ", while missing out " << best_errors << " k-mers.\n";
            
            
            
            if (best_highway.length > 0) {
                // record the score
                res.score += best_highway.length;

                // record the lane switching and matching in the cigar string.
                _update_cigar(cigar, current_lane, current_offset, best_highway);

                // move to the end of the best highway.
                current_lane = best_highway.lane;
                current_offset = best_highway.offset + best_highway.length;
                seqan3::debug_stream << "[INFO]\t\tCurrent position is (" << current_lane << ", " << current_offset << ").\n";
            }
            
            // check if we can break out the loop (no more highways or reached the gaol)
            if (best_highway.length == 0 || current_offset == lanes->l2) {
                break;
            }                
        }
        res.CIGAR = cigar.to_string();
        seqan3::debug_stream << "[INFO]\t\tFinal k-matching score: " << res.score << ", CIGAR string: " << res.CIGAR << ".\n";
        return res;

    }


public:
    /**
     * @brief Construct a new greedy aligner object.
     * 
     * @param kmer_length length of k-mers in the de Brujin graph.
     * @param band_width maximum indel allowed.
     * @param max_errors maximum number of errors between each highways.
     */
    greedy_aligner(unsigned int kmer_length, unsigned int band_width, unsigned int max_errors) : aligner() {
        k = kmer_length;
        bw = band_width;
        e = max_errors;
        lanes = new de_bruijn_lanes<READ_LENGTH>(k, bw);
    }

    ~greedy_aligner() {
        delete lanes;
    }

    align_result_t align(seqan3::dna4_vector s1, seqan3::dna4_vector s2) {
        align_result_t res;
        lanes->read(s1, s2);
        lanes->print();
        _find_subpath();
        return res;
    }
};

#endif