#include <iomanip>
#include "beam_cky_parser.h"
#include "beam_cky_parser.cc"
#include "Utils/reader.h"
#include "Utils/common.h"
#include "Utils/codon.h"

// #ifndef CODON_TABLE
// #define CODON_TABLE "./codon_usage_freq_table_human.csv"
// #endif

#ifndef CODING_WHEEL
#define CODING_WHEEL "./coding_wheel.txt"
#endif

using namespace LinearDesign;

template <typename ScoreType, typename IndexType>
bool output_result(const DecoderResult<ScoreType, IndexType>& result, 
        const double duration, const double lambda, const bool is_verbose, 
        const Codon& codon, string& CODON_TABLE, string& utr_seq) {
    string annotated_seq, cds;

    annotated_seq = result.sequence;
    if (utr_seq != "")
        std::transform(utr_seq.begin(), utr_seq.end(),
                       annotated_seq.end() - utr_seq.size(), ::tolower);

    cds = result.sequence.substr(0, result.sequence.size() - utr_seq.size());

    stringstream ss;
    if (is_verbose)
        ss << "Using lambda = " << (lambda / 100.) << "; Using codon frequency table = " << CODON_TABLE << endl;
    ss << "mRNA sequence:  " << annotated_seq << endl;
    ss << "mRNA structure: " << result.structure << endl;
    ss << "mRNA folding free energy: " << std::setprecision(2) << fixed << result.score 
                                        << " kcal/mol; mRNA CAI: " << std::setprecision(3) 
                                        << fixed << codon.calc_cai(cds) << endl;
    if (is_verbose)
        ss << "Runtime: " << duration << " seconds" << endl;
    cout << ss.str() << endl;

    return true;
}

void show_usage() {
    cerr << "echo SEQUENCE | ./lineardesign -l [LAMBDA]" << endl;
    cerr << "OR" << endl;
    cerr << "cat SEQ_FILE_OR_FASTA_FILE | ./lineardesign -l [LAMBDA]" << endl;
}

// Parse basepairing masks
// Format: "start5,end5,start3,end3 start5,end5,start3,end3 ..."
static struct BasepairingMask *
parse_basepairing_masks(char *arg, int &num_masks)
{
    struct BasepairingMask *masks = NULL;
    char *mask_str = strdup(arg);
    char *saveptr = NULL;
    char *token = strtok_r(mask_str, " ", &saveptr);
    int i = 0;

    while (token != NULL) {
        masks = (struct BasepairingMask *)realloc(masks, (i + 1) * sizeof(struct BasepairingMask));
        if (masks == NULL) {
            free(mask_str);
            return NULL;
        }

        if (sscanf(token, "%d,%d,%d,%d", &masks[i].start5, &masks[i].end5,
                   &masks[i].start3, &masks[i].end3) != 4) {
            free(mask_str);
            free(masks);
            return NULL;
        }

        i++;
        token = strtok_r(NULL, " ", &saveptr);
    }

    free(mask_str);
    num_masks = i;

    return masks;
}

int main(int argc, char** argv) {

    // default args
    double lambda = 0.0f;
    bool is_verbose = false;
    string CODON_TABLE = "./codon_usage_freq_table_human.csv";
    struct PLDOptions pld_options;

    memset(&pld_options, 0, sizeof(pld_options));

    if (getenv("_PLD_REQUEST") != NULL)
        cout << "PLD-capability:" << PLD_CAPABILITY << endl;

    // parse args
    if (argc != 7) {
        show_usage();
        return 1;
    }else{
        lambda = atof(argv[1]);
        is_verbose = atoi(argv[2]) == 1;
        if (string(argv[3]) != ""){
            CODON_TABLE = argv[3];
        }
        pld_options.basepairing_masks = parse_basepairing_masks(
            argv[4], pld_options.num_basepairing_masks);
        srandom((unsigned int)atoi(argv[5]));
        pld_options.random_rejection_threshold =
        pld_options.random_rejection_effective_threshold = RAND_MAX * atof(argv[6]);
    } 
    lambda *= 100.;
    
    // load codon table and coding wheel
    Codon codon(CODON_TABLE);
    std::unordered_map<string, Lattice<IndexType>> aa_graphs_with_ln_weights;
    std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    prepare_codon_unit_lattice<IndexType>(CODING_WHEEL, codon, aa_graphs_with_ln_weights, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon, lambda);

    // start design
    do {
        string aa_seq, aa_tri_seq, utr_seq;

        // convert to uppercase
        transform(aa_seq.begin(), aa_seq.end(), aa_seq.begin(), ::toupper);

        DFA<IndexType> dfa = get_dfa<IndexType>(aa_seq, utr_seq);

        if (is_verbose)
            cout << "Input protein: " << aa_seq << endl;
        if (!ReaderTraits<Fasta>::cvt_to_seq(aa_seq, aa_tri_seq)) 
            continue;

        // init parser
        BeamCKYParser<ScoreType, IndexType> parser(lambda, is_verbose, &pld_options);

        auto protein = util::split(aa_tri_seq, ' ');
        // parse
        auto system_start = chrono::system_clock::now();
        auto result = parser.parse(dfa, codon, aa_seq, protein, aa_best_path_in_a_whole_codon, best_path_in_one_codon_unit);
        auto system_diff = chrono::system_clock::now() - system_start;
        auto system_duration = chrono::duration<double>(system_diff).count();  

        // output
        output_result(result, system_duration, lambda, is_verbose, codon, CODON_TABLE, utr_seq);

#ifdef FINAL_CHECK
        if (codon.cvt_rna_seq_to_aa_seq(result.sequence) != aa_seq) {
            std::cerr << "Final Check Failed:" << std::endl;
            std::cerr << codon.cvt_rna_seq_to_aa_seq(result.sequence) << std::endl;
            std::cerr << aa_seq << std::endl;
            assert(false);
        }
#endif
    } while (0);

    if (pld_options.basepairing_masks != NULL) {
        free(pld_options.basepairing_masks);
    }

    return 0;
}
