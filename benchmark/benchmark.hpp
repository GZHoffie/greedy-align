#ifndef GREEDY_ALIGN_BENCHMARK_HPP
#define GREEDY_ALIGN_BENCHMARK_HPP

#include "../utils.hpp"

#include "../aligner/greedy_aligner.hpp"
#include "../aligner/k_matching.hpp"
#include "../aligner/seqan_aligner.hpp"

#include "../tools/pairwise_dataset.hpp"

template <unsigned int READ_LENGTH>
class aligner_benchmark {
private:
    // aligners
    greedy_aligner<READ_LENGTH>* greedy_align;
    exact_k_matching_aligner<READ_LENGTH>* k_matching_align;
    seqan_aligner* seqan_align;

    std::vector<int> greedy_align_scores;
    std::vector<int> k_matching_align_scores;
    std::vector<int> seqan_align_scores;

    // dataset generation
    read_text_pair_simulator* simulator;
    std::filesystem::path output_data_path;
    std::string indicator_str;

    std::string _get_dataset_filename(float substitution_rate, float indel_rate) {
        return indicator_str + "_" + std::to_string(substitution_rate) + "_" + std::to_string(indel_rate);
    }

    std::string _generate_dataset(int read_len, float substitution_rate, float indel_rate, unsigned int num_reads) {
        simulator->set(read_len, substitution_rate, indel_rate, indel_rate);
        std::string filename = _get_dataset_filename(substitution_rate, indel_rate);
        simulator->generate_dataset_file(output_data_path,
                                         filename,
                                         num_reads);
        return filename;
    }

public:
    aligner_benchmark(unsigned int kmer_length, unsigned int band_width, unsigned int max_errors) {
        // initialize the aligners
        greedy_align = new greedy_aligner<READ_LENGTH>(kmer_length, band_width, max_errors);
        k_matching_align = new exact_k_matching_aligner<READ_LENGTH>((uint8_t)kmer_length);
        seqan_align = new seqan_aligner;

        simulator = new read_text_pair_simulator(kmer_length);
    }

    ~aligner_benchmark() {
        delete greedy_align;
        delete k_matching_align;
        delete seqan_align;
        delete simulator;
    }

    /**
     * @brief Initialize the simulator.
     * 
     * @param fasta_file_name name of the reference genome file.
     * @param indicator the name of the generated
     */
    void initialize(std::filesystem::path const & fasta_file_name,
                    std::filesystem::path const & output_dataset_path, 
                    std::string indicator) {
        simulator->read(fasta_file_name);
        output_data_path = output_dataset_path;
        indicator_str = indicator;
    }

    /**
     * @brief Test the aligner on the dataset.
     * 
     * @param align_alg pointer to the aligner.
     * @param indicator the string indicating name of the dataset.
     * @param res_vector vector storing the scores of alignments.
     */
    void test_aligner(aligner* align_alg, std::string indicator, std::vector<int>& res_vector) {
        // clear everything in the vector
        res_vector.clear();

        // process the dataset
        seqan3::sequence_file_input<_dna4_traits> fin_read{output_data_path / (indicator + "_read.fasta")};
        seqan3::sequence_file_input<_dna4_traits> fin_text{output_data_path / (indicator + "_text.fasta")};

        // Initialize timers
        Timer timer;
        timer.tick();

        // process reads in the dataset
        for (auto && [text, read] : seqan3::views::zip(fin_text, fin_read)) {
            auto res = align_alg->align(text.sequence(), read.sequence());
            res_vector.push_back(res.score);
        }
        timer.tock();

        // benchmark the time usage
        seqan3::debug_stream << "[BENCHMARK]\tTime used by " << align_alg->name << ": " << timer.elapsed_seconds() << " s (" 
                             << (float)timer.elapsed_seconds() / res_vector.size() * 1000 * 1000 << " μs per pairwise alignment).\n";

        //seqan3::debug_stream << res_vector << "\n";
        
        // benchmark correctness for greedy/k-matching aligner
        if (align_alg->name == "Greedy Aligner" || align_alg->name == "K-Matching Aligner") {
            // read the ground truth file
            std::ifstream is(output_data_path / (indicator + ".ground_truth"));

            std::string cigar;
            int k_matching;

            unsigned int index = 0;

            // benchmark results
            unsigned int correct = 0;
            unsigned int absolute_error = 0;
            int max_error = 0;

            while (is >> cigar >> k_matching) {
                if (k_matching == res_vector[index]) {
                    correct++;
                }
                absolute_error += std::abs(k_matching - res_vector[index]);
                max_error = std::max(max_error, std::abs(k_matching - res_vector[index]));
                index++;
            }

            // output correctness
            seqan3::debug_stream << "[BENCHMARK]\tCorrect K-matching output by " << align_alg->name << ": " << correct << " (" 
                                 << (float)correct / res_vector.size() * 100 << " \% correct rate).\n";
            unsigned int num_errors = res_vector.size() - correct;
            if (num_errors == 0) num_errors = 1;
            seqan3::debug_stream << "[BENCHMARK]\tMean absolute error by " << align_alg->name << ": "
                                 << (float)absolute_error / res_vector.size() << " per alignment. ("
                                 << (float)absolute_error / num_errors << " per inaccurate alignment)\n";
            seqan3::debug_stream << "[BENCHMARK]\tMax absolute error by " << align_alg->name << ": "
                                 << max_error << ".\n";
        }
    }

    void test(int read_len,
              std::vector<float> substitution_rates,
              std::vector<float> indel_rates,
              unsigned int num_reads
              ) {
        // initialize the result storing vectors
        greedy_align_scores.reserve(num_reads);
        k_matching_align_scores.reserve(num_reads);
        seqan_align_scores.reserve(num_reads);


        for (unsigned int i = 0; i < substitution_rates.size(); i++) {
            // generate dataset and output in files.
            auto indicator = _generate_dataset(read_len, substitution_rates[i], indel_rates[i], num_reads);

            // test different aligners with the datasets
            test_aligner(greedy_align, indicator, greedy_align_scores);
            test_aligner(seqan_align, indicator, seqan_align_scores);
            test_aligner(k_matching_align, indicator, k_matching_align_scores);
        }
    }

};

#endif