#ifndef GREEDY_ALIGN_BENCHMARK_TRACE_RECON_HPP
#define GREEDY_ALIGN_BENCHMARK_TRACE_RECON_HPP

#include "../utils.hpp"
#include "../ensembler/pog_ensembler.hpp"

template <unsigned int READ_LENGTH>
class trace_reconstruction_benchmark {
private:
    // ensemblers
    partial_order_graph_ensembler<READ_LENGTH>* pog_ensembler;

    // dataset directory
    std::filesystem::path clusters_path;
    std::filesystem::path centers_path;

    std::vector<seqan3::dna4_vector> _string_vector_to_dna4(const std::vector<std::string>& sequences) {
        std::vector<seqan3::dna4_vector> res;
        for (auto& s : sequences) {
            seqan3::dna4_vector dna_sequence;
            for (auto &nt : s) {
                seqan3::dna4 dna = seqan3::assign_char_to(nt, seqan3::dna4{});
                dna_sequence.push_back(dna);
            }
            res.push_back(dna_sequence);
            //seqan3::debug_stream << res.size() << " " << sequences.size() << "\n";
        }
        return res;
    }

    std::string _dna4_to_string(const seqan3::dna4_vector& sequence) {
        std::string res;
        for (auto & nt : sequence) {
            res += nt.to_char();
        }
        return res;
    }

public:
    trace_reconstruction_benchmark(std::filesystem::path clusters_file, std::filesystem::path centers_file, bool debug = false) {
        pog_ensembler = new partial_order_graph_ensembler<READ_LENGTH>(6, 4, 7, 5, 110, debug);
        clusters_path = clusters_file;
        centers_path = centers_file;
    }

    ~trace_reconstruction_benchmark() {
        delete pog_ensembler;
    }

    void benchmark() {
        std::ifstream clusters_stream(clusters_path);
        std::ifstream centers_stream(centers_path);

        std::string center;
        
        // ignore the first line in the clusters file
        clusters_stream >> center;

        std::vector<std::string> cluster;

        unsigned int total = 0;
        unsigned int correct = 0;

        while (centers_stream >> center) {
            //seqan3::debug_stream << center;
            cluster.clear();

            seqan3::debug_stream << total << " " << center << "\n";
            
            // read the cluster of sequences
            std::string temp;
            while (clusters_stream >> temp) {
                //seqan3::debug_stream << temp << "\n";
                if (temp[0] == '=') break;
                cluster.push_back(temp);
            }

            // turn the cluster to seqan3::dna4_vector's.
            auto dna4_cluster = _string_vector_to_dna4(cluster);

            // ensembling
            auto pog_result = pog_ensembler->ensemble(dna4_cluster);
            if (_dna4_to_string(pog_result) == center) {
                correct++;
            } else {
                seqan3::debug_stream << pog_result << "\n";
                seqan3::debug_stream << center << "\n";
            }
            
            total++;
        }
        seqan3::debug_stream << "Total: " << total << ", correct: " << correct << "\n";
    }
};


#endif