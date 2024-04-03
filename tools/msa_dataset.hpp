#ifndef GREEDY_ALIGN_MSA_DATASET_HPP
#define GREEDY_ALIGN_MSA_DATASET_HPP

#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/alphabet/all.hpp>

#include <cstdlib>
#include <ctime>
#include <random>
#include <tuple>

#include "../utils.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_cigar_operation;

class read_text_msa_simulator {
    /**
     * @brief A tool to randomly generate short reads based on a reference genome.
     *        It allows setting of substitution/indel error rates and random generation
     *        of errors.
     */
private:
    // variables storing the sequence information
    std::vector<std::vector<seqan3::dna4>> sequences;
    int read_length;
    int reference_length;

    // Error generation
    std::vector<seqan3::dna4> neucleotides{ 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4 };
    std::poisson_distribution<int>* read_substitution_dist; 
    std::poisson_distribution<int>* read_indel_dist;
    float SNP_dist;
    float indel_dist;

    // Random number generator
    std::mt19937* gen;

    // Checking the k-matching between the two strings.
    unsigned int k;

    // sample a long sequence from the reference genome
    seqan3::dna4_vector reference_sequence;
    seqan3::dna4_vector consensus_ground_truth;

    // Error generating functions
    void _add_substitution(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        unsigned int index = rand() % sequence.size();
        seqan3::dna4 orig_nt = sequence[index];
        seqan3::dna4 new_nt = neucleotides[rand() % 4];
        while (new_nt == orig_nt) {
            new_nt = neucleotides[rand() % 4];
        }
        sequence[index] = new_nt;
        cigar.replace(index, 'X'_cigar_operation);
    }

    void _add_insertion(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        unsigned int index = rand() % sequence.size();
        sequence.insert(sequence.begin() + index, neucleotides[rand() % 4]);
        cigar.insert(index, 'I'_cigar_operation);
    }

    void _add_deletion(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        unsigned int index = rand() % sequence.size();
        sequence.erase(sequence.begin() + index);
        cigar.replace(index, 'D'_cigar_operation);
    }

public:
    read_text_msa_simulator(int reference_len, int read_len, float SNP_rate = 0, float indel_rate = 0, float read_substitution_rate = 0, 
                            float read_indel_rate = 0, unsigned int seed_length = 3) {
        // set a random seed
        srand(time(NULL));
        
        // Set indexing-related numbers
        reference_length = reference_len;
        read_length = read_len;
        
        // Set the error rates
        read_substitution_dist = new std::poisson_distribution<int>(read_substitution_rate * read_length);
        read_indel_dist = new std::poisson_distribution<int>(read_indel_rate / 2 * read_length);

        SNP_dist = SNP_rate;
        indel_dist = indel_rate;

        // set value of k in k-matching
        k = seed_length;
    }

    /**
     * @brief Constructor, but determine the parameters later.
     * 
     */
    read_text_msa_simulator(unsigned int seed_length = 3) : read_text_msa_simulator(0, 0, 0, 0, seed_length) {}

    ~read_text_msa_simulator() {
        delete read_substitution_dist, read_indel_dist;
    }

    void add_errors(std::vector<seqan3::dna4>& sequence, CIGAR& cigar, 
                    int substitutions, int deletions, int insertions) {
        /**
         * @brief Add errors to a short read sequence.
         * @param sequence the short read sequence to be modified.
         * @param cigar the ground truth cigar string.
         * @param substitutions number of substitutions to make.
         * @param deletions number of deletions to make.
         * @param insertions number of insertions to make.
         */
        for (int i = 0; i < deletions; i++) _add_deletion(sequence, cigar);
        for (int i = 0; i < insertions; i++) _add_insertion(sequence, cigar);
        for (int i = 0; i < substitutions; i++) _add_substitution(sequence, cigar);
    }

    void simulate_read_errors(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        /**
         * @brief Do simulation and generate the number of substitution, deletion and insertion randomly
         *        from poisson distributions with lamda = substitution_rate, deletion_rate and 
         *        insertion_rate, respectively.
         * @note the number of errors generated might always be 0 if error rates are too small. Might consider
         *       setting them higher than actual error rate.
         */
        std::random_device rd;
        std::mt19937 gen(rd());
        add_errors(sequence, cigar,
                   (*read_substitution_dist)(gen),
                   (*read_indel_dist)(gen),
                   (*read_indel_dist)(gen));
    }

    void read(std::filesystem::path const & fasta_file_name) {
        /**
         * @brief Read the fasta file and store the sequences for each bucket. The reference
         *        genome is then used for short read generation.
         * @param fasta_file_name the name of the sequence file to be read.
         */

        seqan3::sequence_file_input<_dna4_traits> fin{fasta_file_name};
 
        for (auto & record : fin) {
            sequences.push_back(record.sequence());
            // TODO: also consider the quality?
        }
    }


    CIGAR sample_reference() {
        /**
         * @brief sample a part of length `reference_length` from the reference genome, and randomly add 
         *        SNPs and indels and store in `consensus_ground_truth`, from where the reads will be sampled.
         * @returns a CIGAR string indicating places of indels and SNPs.
         */
        int sequence = rand() % sequences.size();
        std::vector<seqan3::dna4>& current_sequence = sequences[sequence];
        int size = current_sequence.size();
        int start = 0;
        if (current_sequence.size() > reference_length + 1) {
            start = rand() % (current_sequence.size() - reference_length - 1);
        }
        int end = start + reference_length;
        if (end > current_sequence.size()) {
            end = current_sequence.size();
        }
        std::vector<seqan3::dna4> sample_sequence(current_sequence.begin() + start, current_sequence.begin() + end);

        // copy the original sequence
        reference_sequence = sample_sequence;
        consensus_ground_truth = sample_sequence;

        // add some errors to be the consensus
        int length = reference_sequence.size();
        CIGAR cigar(length, '='_cigar_operation);
        add_errors(consensus_ground_truth, cigar, (int)(SNP_dist * length), (int)(indel_dist / 2 * length), (int)(indel_dist / 2 * length));

        return cigar;
    }

    std::tuple<int, std::vector<seqan3::dna4>, CIGAR> sample() {
        /**
         * @brief Take a sample short read from the reference genome and add errors to it.
         * @return a tuple storing 3 values: <Starting point of the generated read, the original genome,
         *         the starting point of the short read>.
         */
        std::vector<seqan3::dna4>& current_sequence = consensus_ground_truth;
        int size = current_sequence.size();
        int start = 0;
        if (current_sequence.size() > read_length + 1) {
            start = rand() % (current_sequence.size() - read_length - 1);
        }
        int end = start + read_length;
        if (end > current_sequence.size()) {
            end = current_sequence.size();
        }
        std::vector<seqan3::dna4> sample_sequence(current_sequence.begin() + start, current_sequence.begin() + end);

        // copy the original sequence
        std::vector<seqan3::dna4> original_sequence(sample_sequence);
        // assign the original cigar sequence
        CIGAR cigar(sample_sequence.size(), '='_cigar_operation);

        // insert errors
        //simulate_errors(sample_sequence, cigar);
        return std::make_tuple(start, sample_sequence, cigar);
    }

    void generate_dataset_file(std::filesystem::path output_path, std::string indicator, unsigned int size) {
        /**
         * @brief Generate two fasta files containing short reads.
         * @param output_path the directory where we output the fastq file and the answer file.
         * @param indicator a string that indicate the name of output fastq/answer file.
         * @param size the number of short reads to be simulated.
         */
        // Create directory if directory is not created yet.
        // Return if the index files already exist.
        if (!std::filesystem::create_directories(output_path)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified output directory "
                                 << output_path << " is already created." << '\n';
            if (!check_filename_in(output_path, indicator + "_read.fasta") || 
                !check_filename_in(output_path, indicator + "_reference.fasta") ||
                !check_filename_in(output_path, indicator + ".ground_truth")) return;
        }

        std::ofstream text_file(output_path / (indicator + "_reference.fasta"));
        std::ofstream read_file(output_path / (indicator + "_read.fasta"));
        std::ofstream read_ground_truth_file(output_path / (indicator + ".read_ground_truth"));
        std::ofstream mutation_ground_truth_file(output_path / (indicator + ".mutation_ground_truth"));

        // sample the reference and write the ground truth to the corresponding files.
        auto cigar = sample_reference();
        text_file << ">Sampled Reference\n";
        for (auto nt : reference_sequence) {
            text_file << nt.to_char();
        }
        text_file << "\n";
        mutation_ground_truth_file << cigar.to_string() << "\n";


        for (unsigned int i = 0; i < size; i++) {
            auto [start, read, cigar] = sample();
            read_file << ">" << i << " " << start << "\n";

            // write sequence
            for (auto nt : read) {
                read_file << nt.to_char();
            }
            
            read_file << "\n";

            // record the ground truth, including the ground truth cigar string and the k-matching.
            read_ground_truth_file << cigar.to_string() << "\n";
        }

        seqan3::debug_stream << "[INFO]\t\t" << "The generated fasta file is stored in: " 
                             << output_path / (indicator + "_[read/reference].fasta") << ".\n";
    }

};

#endif