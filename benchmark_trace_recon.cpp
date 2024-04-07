
#include "./utils.hpp"

#include "benchmark/benchmark_trace_recon.hpp"

int main() {
    // name/path of the dataset
    std::filesystem::path centers_path = "/home/zhenhao/greedy-align/data/Centers.txt";
    std::filesystem::path clusters_path = "/home/zhenhao/greedy-align/data/Clusters.txt";
    //std::filesystem::path centers_path = "/home/zhenhao/greedy-align/data/Centers_test.txt";
    //std::filesystem::path clusters_path = "/home/zhenhao/greedy-align/data/Clusters_test.txt";

    trace_reconstruction_benchmark<120> bench(clusters_path, centers_path, false);

    bench.benchmark();     
    return 0;
}