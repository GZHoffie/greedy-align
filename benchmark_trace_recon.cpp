
#include "./utils.hpp"

#include "benchmark/benchmark_trace_recon.hpp"

int main() {
    // name/path of the dataset
    //std::filesystem::path centers_path = "/home/zhenhao/greedy-align/data/Centers.txt";
    //std::filesystem::path clusters_path = "/home/zhenhao/greedy-align/data/Clusters.txt";
    //std::filesystem::path centers_path = "/home/zhenhao/greedy-align/data/Centers_test.txt";
    //std::filesystem::path clusters_path = "/home/zhenhao/greedy-align/data/Clusters_test.txt";

    std::filesystem::path centers_path = "/mnt/c/Users/zhenh/trace_recon/our_nanopore_refs.txt";
    std::filesystem::path clusters_path = "/mnt/c/Users/zhenh/trace_recon/our_nanopore_UnderlyingClusters.txt";

    //std::filesystem::path centers_path = "/mnt/c/Users/zhenh/trace_recon/our_illum_refs.txt";
    //std::filesystem::path clusters_path = "/mnt/c/Users/zhenh/trace_recon/our_illum_UnderlyingClusters.txt";

    //std::filesystem::path centers_path = "/mnt/c/Users/zhenh/trace_recon/our_illum_refs.txt";
    //std::filesystem::path clusters_path = "/mnt/c/Users/zhenh/trace_recon/our_illum_UnderlyingClusters.txt";

    std::filesystem::path output_path = "/mnt/c/Users/zhenh/trace_recon/Output_nanopore.txt";

    trace_reconstruction_benchmark<120> bench(clusters_path, centers_path, output_path, false);

    bench.benchmark();     
    return 0;
}