#include "aligner/greedy_aligner.hpp"
#include "aligner/k_matching.hpp"
#include "tools/generate_dataset.hpp"
#include "./utils.hpp"

int main() {
    // parameters
    int k = 3;
    int errors = 5;
    int bw = 2;

    // name/path of the dataset
    std::filesystem::path genome_path = "/mnt/d/genome/GRCh38_latest_genomic.fasta";
    std::filesystem::path dataset_path = "/mnt/d/genome/align";
    std::string indicator = "GRCh38_150";

    // settings of the generated dataset
    int read_length = 150;
    int read_number = 10;

    // initialize aligners
    exact_k_matching_aligner<160> k_matching_align(k);
    greedy_aligner<160> greedy_align(k, bw, errors);


    // generate dataset
    bool generate_dataset = false;
    if (generate_dataset) {
        read_text_pair_simulator simulator(read_length, 0.002, 0.00025, 0.00025);
        simulator.read(genome_path);
        simulator.generate_dataset_file(dataset_path, indicator, read_number);
    }

    // process the dataset
    seqan3::sequence_file_input<_dna4_traits> fin_read{dataset_path / (indicator + "_read.fasta")};
    seqan3::sequence_file_input<_dna4_traits> fin_text{dataset_path / (indicator + "_text.fasta")};

    int correct = 0;
    int index = 0, test_num = 1000;
 
    for (auto && [text, read] : seqan3::views::zip(fin_text, fin_read)) {
        auto res_greedy = greedy_align.align(text.sequence(), read.sequence());
        auto res_exact = k_matching_align.align(text.sequence(), read.sequence());
        
        if (res_greedy.score == res_exact.score) {
            correct++;
        } else {
            seqan3::debug_stream << res_greedy.score << " | " << res_exact.score << "\n";
            seqan3::debug_stream << text.sequence() << "\n";
            seqan3::debug_stream << read.sequence() << "\n";
        }
        index++;
        if (index >= test_num) break;
    }

    seqan3::debug_stream << "Correct: " << correct << "\n";
    return 0;
}