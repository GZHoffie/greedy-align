#include "aligner.hpp"

#include <seqan3/search/kmer_index/shape.hpp>
#include <bitset>

template<unsigned int READ_LENGTH>
class de_bruijn_lanes {
private:
    std::vector<std::bitset<READ_LENGTH>> lanes;
    unsigned int bw;

public:
    /**
     * @brief Construct a new de_bruijn_lanes object.
     * 
     * @param band_width Maximum number of indels allowed.
     */
    de_bruijn_lanes(int band_width) {
        bw = band_width;
        for (int i = -band_width; i <= band_width; i++) {
            std::bitset<READ_LENGTH> lane;
            lanes.push_back(lane);
        }
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
    std::bitset<READ_LENGTH>& at(int i) {
        return lanes[i + bw];
    }

    /**
     * @brief Print out the lanes in order, for debugging.
     * 
     */
    void print() {
        seqan3::debug_stream << "[INFO]\t\t" << "Current de Bruijn lane layout:\n";
        for (int i = 0; i < lanes.size(); i++) {
            seqan3::debug_stream << "[INFO]\t\t" << "Lane " << i- (int)bw << "\t| " << lanes.at(i) << "\n";
        }
    }
};

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

template<unsigned int READ_LENGTH>
class greedy_aligner : public aligner {
private:
    unsigned int k;
    int bw;
    de_bruijn_lanes<READ_LENGTH>* lanes;


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
        lanes->at(shift_amount) = v1.shifted_xor(v2, shift_amount);
        
    }

    /**
     * @brief Find highways (common subpaths) in the de Bruijn graphs, which are represented by
     *        (consecutive) zeros in the lanes.
     * 
     * @param s1, s2 the DNA sequences represented using seqan3::dna4 objects.
     */
    void _build_de_bruijn_path(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) {
        // build read_bit_vectors objects base on the reads
        read_bit_vectors<READ_LENGTH> v1(s1);
        read_bit_vectors<READ_LENGTH> v2(s2);

        // perform shifting and bitwise xor operation, and store in `lanes`.
        for (int l = -bw; l <= bw; l++) {
            _shifted_xor(v1, v2, l);
            auto original_lane = lanes->at(l);
            for (int i = 1; i < k; i++) {
                lanes->at(l) |= (original_lane >> i);
            }
        }
    }


public:
    /**
     * @brief Construct a new greedy aligner object.
     * 
     * @param kmer_length length of k-mers in the de Brujin graph.
     * @param band_width maximum indel allowed.
     */
    greedy_aligner(unsigned int kmer_length, unsigned int band_width) : aligner() {
        k = kmer_length;
        bw = band_width;
        lanes = new de_bruijn_lanes<READ_LENGTH>(bw);
    }

    ~greedy_aligner() {
        delete lanes;
    }

    align_result_t align(seqan3::dna4_vector s1, seqan3::dna4_vector s2) {
        align_result_t res;
        _build_de_bruijn_path(s1, s2);
        lanes->print();
        return res;
    }
};