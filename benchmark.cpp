#include "aligner/greedy_aligner.hpp"
#include "aligner/parasail_aligner.hpp"
#include "aligner/seqan_aligner.hpp"
#include "aligner/k_matching.hpp"
#include "tools/generate_dataset.hpp"
#include "./utils.hpp"

int main() {
    // parameters
    int k = 5;
    int errors = 5;
    int bw = 4;

    // name/path of the dataset
    std::filesystem::path genome_path = "/mnt/d/genome/GRCh38_latest_genomic.fasta";
    std::filesystem::path dataset_path = "/mnt/d/genome/align";
    std::string indicator = "GRCh38_150";

    // settings of the generated dataset
    int read_length = 150;
    int read_number = 1000000;

    // initialize aligners
    exact_k_matching_aligner<160> k_matching_align(k);
    seqan_aligner greedy_align;


    // generate dataset
    bool generate_dataset = false;
    if (generate_dataset) {
        read_text_pair_simulator simulator(read_length, 0.1, 0.05, 0.05);
        simulator.read(genome_path);
        simulator.generate_dataset_file(dataset_path, indicator, read_number);
    }

    // process the dataset
    seqan3::sequence_file_input<_dna4_traits> fin_read{dataset_path / (indicator + "_read.fasta")};
    seqan3::sequence_file_input<_dna4_traits> fin_text{dataset_path / (indicator + "_text.fasta")};

    int correct = 0;
    int index = 0, test_num = 100000;

    // Initialize timers
    Timer greedy_timer;
    Timer exact_matching_timer;
 
    greedy_timer.tick();
    for (auto && [text, read] : seqan3::views::zip(fin_text, fin_read)) {
        auto res_greedy = greedy_align.align(text.sequence(), read.sequence());
        /*
        auto res_exact = k_matching_align.align(text.sequence(), read.sequence());
        if (res_greedy.score == res_exact.score) {
            correct++;
        } else {
            //seqan3::debug_stream << res_greedy.score << " | " << res_exact.score << "\n";
            //seqan3::debug_stream << text.sequence() << "\n";
            //seqan3::debug_stream << read.sequence() << "\n";
        }
        index++;
        if (index >= test_num) break;
        */
        index++;
    }
    greedy_timer.tock();

    seqan3::debug_stream << "Correct: " << correct << "\n";
    seqan3::debug_stream << "[BENCHMARK]\tTime used by greedy aligner: " << greedy_timer.elapsed_seconds() << " s (" 
                         << (float)greedy_timer.elapsed_seconds() / index * 1000 * 1000 << " μs per pairwise alignment).\n";                     
    return 0;
}