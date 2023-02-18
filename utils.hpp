#ifndef GREEDY_ALIGN_UTILS_HPP
#define GREEDY_ALIGN_UTILS_HPP


#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

#include <chrono>
#include <thread>
#include <assert.h>

/**
 * @brief Class for benchmark of runtime.
 *        Adopted from https://coliru.stacked-crooked.com/a/508ec779dea5a28d.
 */
template <class DT = std::chrono::milliseconds,
          class ClockT = std::chrono::steady_clock>
class Timer
{
    using timep_t = decltype(ClockT::now());
    
    // variable storing start and end time.
    timep_t _start = ClockT::now();
    timep_t _end = {};

    // total duration between each (_start, _end) pair
    DT _elapsed = DT::zero();

public:
    auto duration() const { 
        // Use gsl_Expects if your project supports it.
        assert(_end != timep_t{} && "Timer must toc before reading the time"); 
        return _elapsed; 
    }

    void tick() { 
        _end = timep_t{};
        _start = ClockT::now(); 
    }
    
    void tock() {
        _end = ClockT::now();
        _elapsed += std::chrono::duration_cast<DT>(_end - _start);
    }
    
    float elapsed_seconds() const {
        return ((float) _elapsed.count()) / 1000;
    }
};


class CIGAR {
    /**
     * @brief Util class for the cigar string, with support of inserting cigar operation.
     */
private:
    std::vector<seqan3::cigar> CIGAR_vector;

public:
    CIGAR() = default;

    explicit CIGAR(std::vector<seqan3::cigar> CIGAR_init) {
        /**
         * @brief Initialize the class with a vector of seqan3::cigar objects.
         * 
         */
        CIGAR_vector = CIGAR_init;
    }

    /**
     * @brief Initialize the class with `size` number of `op`s.
     * 
     * @param size The number of operations in the CIGAR string.
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc.
     */
    CIGAR(unsigned int size, seqan3::cigar::operation op) {
        /**
         * @brief Initialize the class with `size` number of `op`s.
         * 
         */
        for (int i = 0; i < size; i++) {
            seqan3::cigar operation{1, op};
            CIGAR_vector.push_back(operation);
        }
    }

    /**
     * @brief Insert a cigar operation at a certain index.
     * 
     * @param size The number of operations in the CIGAR string.
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc.
     */
    void insert(unsigned int index, seqan3::cigar::operation op) {
        seqan3::cigar operation{1, op};
        CIGAR_vector.insert(CIGAR_vector.begin() + index, operation);
    }

    /**
     * @brief Insert a cigar operation at the back of all operations.
     * 
     * @param size The number of operations in the CIGAR string.
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc.
     */
    void push_back(unsigned int size, seqan3::cigar::operation op) {
        seqan3::cigar operation{size, op};
        CIGAR_vector.push_back(operation);
    }

    /**
     * @brief Replace the operation at `index` to `op`.
     * 
     * @param index the index at which the operation is replaced
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc. 
     */
    void replace(unsigned int index, seqan3::cigar::operation op) {
        seqan3::cigar operation{1, op};
        CIGAR_vector[index] = operation;
    }

    /**
     * @brief Convert the vector to a CIGAR string for output.
     * 
     * @return std::string the CIGAR string that captures the whole vector.
     */
    std::string to_string() {
        std::string res;
        unsigned int size = 0;
        seqan3::cigar::operation last_op, curr_op;
        if (!CIGAR_vector.empty()) {
            last_op = get<1>(CIGAR_vector[0]);
            for (auto cigar : CIGAR_vector) {
                auto &[curr_size, curr_op] = cigar;
                if (last_op == curr_op) {
                    size += curr_size;
                } else if (size != 0) {
                    res += std::to_string(size) + last_op.to_char();
                    size = curr_size;
                    last_op = curr_op;
                }
            }
            if (size != 0) {
                res += std::to_string(size) + last_op.to_char();
            }
        }
        return res;
    }
};

#endif