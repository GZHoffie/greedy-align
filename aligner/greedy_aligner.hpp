#include "aligner.hpp"

#include <seqan3/search/kmer_index/shape.hpp>

template<unsigned int READ_LENGTH>
class de_bruijn_lanes {
private:
    std::vector<std::bitset<READ_LENGTH>> lanes;

public:
    /**
     * @brief Construct a new de_bruijn_lanes object.
     * 
     * @param band_width Maximum number of indels allowed.
     */
    de_bruijn_lanes(unsigned int band_width) {
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
        return lanes[i + band_width];
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
            b1[i] = s.to_rank() & 1;
            b2[i] = s.to_rank() & 2;
        }
        length = s.size();
    }

    std::bitset<READ_LENGTH> shifted_xor(const read_bit_vectors& that, int shift_amount) {
        std::bitset<READ_LENGTH> res;
        // void out all meaningless bits
        res = ~res << (READ_LENGTH - that.length);
        if (shift_amount > 0) {
            // void out the shifted bits
            res = ~(res >> (READ_LENGTH - that.length) << shift_amount >> shift_amount);
            // perform bitwise xor operation
            res |= ((b1 << shift_amount) ^ that.b1) | ((b2 << shift_amount) ^ that.b2);
        } else {
            // same as above, but in opposite direction
            res = ~(res << shift_amount >> (shift_amount + READ_LENGTH - that.length));
            res |= ((b1 >> -shift_amount) ^ that.b1) | ((b2 >> -shift_amount) ^ that.b2);
        }
        return res;
    }

};

template<unsigned int READ_LENGTH>
class greedy_aligner : public aligner {
private:
    unsigned int k;
    unsigned int bw;
    de_bruijn_lanes<READ_LENGTH>* lanes;


    /**
     * @brief Check the hamming distance between (v1 << shift_amount) ^ v2, and store
     *        the results in `lanes`.
     * 
     * @param v1, v2 the two sequences, already converted to bit vectors.
     * @param shift_amount The amount of shifting, which should be in [-band_width, band_width].
     */
    void _shifted_xor(const read_bit_vectors& v1, const read_bit_vectors& v2, int shift_amount) {
        lanes->at(shift_amount) = v1.shifted_xor(v2, shift_amount);
    }

    void _build_de_bruijn_path(const seqan3::dna4_vector& s1, const seqan3::dna4_vector& s2) {
        // build read_bit_vectors objects base on the reads
        read_bit_vectors v1(s1);
        read_bit_vectors v2(s2);

        // perform shifting and bitwise xor operation, and store in `lanes`.
        for (int l = -bw; l <= bw; l++) {
            _shifted_xor(v1, v2, l);
            for (int i = 0; i < k; i++) {
                lanes->at(l) |= (lanes->at(l) << i);
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
    }

    align_result_t align(seqan3::dna4_vector s1, seqan3::dna4_vector s2) {
        align_result_t res;
        _build_de_bruijn_path(s1, s2);
        return res;
    }



};