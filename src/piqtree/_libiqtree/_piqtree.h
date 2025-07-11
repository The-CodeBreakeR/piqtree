#ifndef _PIQTREE_H
#define _PIQTREE_H

#include <pybind11/stl.h>
#include <cstddef>

using namespace std;

#ifdef _MSC_VER
#pragma pack(push, 1)
#else
#pragma pack(1)
#endif

typedef struct {
  const char** strings;
  size_t length;
} StringArray;

typedef struct {
  double* doubles;
  size_t length;
} DoubleArray;

typedef struct {
  int value;
  char* errorStr;
} IntegerResult;

typedef struct {
  char* value;
  char* errorStr;
} StringResult;

typedef struct {
  double* value;
  size_t length;
  char* errorStr;
} DoubleArrayResult;

#ifdef _MSC_VER
#pragma pack(pop)
#else
#pragma pack()
#endif

/*
 * Calculates the robinson fould distance between two trees
 */
extern "C" IntegerResult robinson_fould(const char* ctree1, const char* ctree2);

/*
 * Generates a set of random phylogenetic trees
 * tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED",
 * "BIRTH_DEATH", "STAR_TREE" output: a newick tree (in string format)
 */
extern "C" StringResult random_tree(int num_taxa,
                                    const char* tree_gen_mode,
                                    int num_trees,
                                    int rand_seed = 0);

/*
 * Perform phylogenetic analysis on the input alignment
 * With estimation of the best topology
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection
 * of the optimal number of cpu threads output: results in YAML format with the
 * tree and the details of parameters
 */
extern "C" StringResult build_tree(StringArray& names,
                                   StringArray& seqs,
                                   const char* model,
                                   int rand_seed = 0,
                                   int bootstrap_rep = 0,
                                   int num_thres = 1);

/*
 * Perform phylogenetic analysis on the input alignment
 * With restriction to the input toplogy
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection
 * of the optimal number of cpu threads output: results in YAML format with the
 * details of parameters
 */
extern "C" StringResult fit_tree(StringArray& names,
                                 StringArray& seqs,
                                 const char* model,
                                 const char* intree,
                                 int rand_seed = 0,
                                 int num_thres = 1);

/*
 * Perform phylogenetic analysis with ModelFinder
 * on the input alignment (in string format)
 * model_set -- a set of models to consider
 * freq_set -- a set of frequency types
 * rate_set -- a set of RHAS models
 * rand_seed -- random seed, if 0, then will generate a new random seed
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection
 * of the optimal number of cpu threads output: modelfinder results in YAML
 * format
 */
extern "C" StringResult modelfinder(StringArray& names,
                                    StringArray& seqs,
                                    int rand_seed = 0,
                                    const char* model_set = "",
                                    const char* freq_set = "",
                                    const char* rate_set = "",
                                    int num_thres = 1);

/*
 * Build pairwise JC distance matrix
 * output: set of distances
 * (n * i + j)-th element of the list represents the distance between i-th and
 * j-th sequence, where n is the number of sequences num_thres -- number of cpu
 * threads to be used, default: 1; 0 - use all available cpu threads on the
 * machine
 */
extern "C" DoubleArrayResult build_distmatrix(StringArray& names,
                                              StringArray& seqs,
                                              int num_thres = 1);

/*
 * Using Rapid-NJ to build tree from a distance matrix
 * output: a newick tree (in string format)
 */
extern "C" StringResult build_njtree(StringArray& names,
                                     DoubleArray& distances);

/*
 * Compute a consensus tree
 * trees   -- a set of input trees
 * minsup  -- a threshold to build the majority consensus (default: 0.5)
 * output  -- the resulting consensus tree
 */
extern "C" StringResult consensus_tree(StringArray& trees, double minsup = 0.5);

/*
 * verion number
 */
extern "C" StringResult version();

/*
 * Execute AliSim Simulation
 * output: results in YAML format that contains the simulated alignment and the content of the log file
 * trees -- array of NEWICK tree strings (multiple trees)
 * subst_model -- the substitution model name
 * seed -- the random seed
 * partition_info -- partition information
 * partition_type -- partition type is either ‘equal’, ‘proportion’, or ‘unlinked’
 * seq_length -- the length of sequences
 * insertion_rate -- the insertion rate
 * deletion_rate -- the deletion rate
 * root_seq -- the root sequence
 * num_threads -- the number of threads
 * insertion_size_distribution -- the insertion size distribution
 * deletion_size_distribution -- the deletion size distribution
 */
extern "C" StringResult simulate_alignment(StringArray& trees,
                                           const char* subst_model,
                                           int seed,
                                           StringArray& partition_info,
                                           const char* partition_type = "",
                                           int seq_length = 1000,
                                           double insertion_rate = 0,
                                           double deletion_rate = 0,
                                           const char* root_seq = "",
                                           int num_threads = 1,
                                           const char* insertion_size_distribution = "",
                                           const char* deletion_size_distribution = "");


#endif /* _PIQTREE_H */
