
#include "./utils.hpp"

#include "benchmark/benchmark_pairwise.hpp"

int main() {
    // parameters
    int k = 3;
    int errors = 5;
    int bw = 4;

    // name/path of the dataset
    std::filesystem::path genome_path = "/mnt/d/genome/GRCh38_latest_genomic.fasta";
    std::filesystem::path dataset_path = "/mnt/d/genome/align";
    std::string indicator = "GRCh38";

    // settings of the generated dataset
    int read_length = 150;
    int read_number = 100000;

    // initialize aligners
    exact_k_matching_aligner<160> k_matching_align(k);
    seqan_aligner greedy_align;

    // testing error rate
    std::vector<float> substitution_rates{ 0.002, 0.01, 0.05, 0.3 };
    std::vector<float> indel_rates{ 0.00025, 0.00125, 0.00625, 0.0375 };

    // set up benchmark
    aligner_benchmark<160> benchmark(k, bw, errors);   
    benchmark.initialize(genome_path, dataset_path, indicator);
    benchmark.test(read_length, substitution_rates, indel_rates, read_number);              
    return 0;
}