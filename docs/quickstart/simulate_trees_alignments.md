# Simulate coalescent trees with msprime and an alignment with AliSim

A set of coalescent trees can be simulated using `sim_ancestry()` from the [msprime](https://github.com/tskit-dev/msprime). Then, an alignment can be simulated using `simulate_alignment()` from the `piqtree` package.

## Usage

### Example for Single Population

The following pipeline simulates a set of trees independently under the coalescent model with a given population size and recombination rate using `msprime`, then rescales the branch lengths of the trees to substitution units, and finally simulates an alignment for the rescaled trees using the `piqtree`'s version of `AliSim` with a given substitution model and random seed.

```python
from piqtree import simulate_alignment
import cogent3
import msprime

# Set simulation parameters
num_taxa = 10
seq_length = 2000
recombination_rate = 1e-8
population_size = 10000
num_threads = 5
seed = 1

# Simulate trees under the coalescent model using msprime
ts = msprime.sim_ancestry(
        samples=num_taxa/2,
        sequence_length=seq_length,
        recombination_rate=recombination_rate,
        population_size=population_size,
        random_seed=seed)

# Take the first tree for simulating an alignment
for t in ts.trees():
    newick_tree = t.as_newick()
    break
tree = cogent3.make_tree(newick_tree)

# Simulate alignment from the rescaled simulated tree using AliSim
aln, log = simulate_alignment(
        trees = [tree],
        model = "JC",
        rand_seed = seed,
        num_threads = num_threads,
        length = seq_length,
        population_size = population_size)

# Print the alignment
print(aln)

# Print the console log
print(log)

```

### Example for Multiple Populations/Species

The following pipeline extends the previous example by simulating a set of gene trees generated under the multi-species coalescent model given a species tree with four species `human`, `chimp`, `gorilla` and `orangutan`. 


```python
import msprime
from piqtree import simulate_alignment
import cogent3

# Set simulation parameters
seq_length = 2000
recombination_rate = 1e-8
population_size = 100_000 # larger population sizes create conditions with more gene tree discordance
generation_time = 20
num_threads = 5
seed = 1
num_genes = 10

species_list = ["human", "chimp", "gorilla", "orangutan"]
num_ind_per_species = [1, 2, 1, 3] # sets varying number of individuals per species
species_tree = "(((human:5.6,chimp:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)"

samples, sample_labels, sample_index = [], {}, 0
for i in range(len(species_list)):
    samples.append(msprime.SampleSet(num_ind_per_species[i], population=species_list[i], time=0))
for i in range(len(species_list)):
    for j in range(2*num_ind_per_species[i]):  # 2 samples per species
        sample_labels[sample_index] = f"{species_list[i]}_{j+1}"
        sample_index += 1

# Create demography from species tree
demog = msprime.Demography.from_species_tree(
    species_tree,
    time_units="myr",
    generation_time=generation_time,
    initial_size=population_size)

gene_trees = []

# Simulate gene trees given the demography under the coalescent model
for i in range(num_genes):
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demog,
        sequence_length=seq_length,
        random_seed=seed+i)
    for t in ts.trees():
        newick = t.as_newick(node_labels=sample_labels)
        gene_trees.append(cogent3.make_tree(newick))
        print(f"Tree {i}:\n{newick}\n")
        break

# Create a partition file
partition_info = '#nexus\n begin sets;\n\t'
for i in range(len(gene_trees)):
    partition_info += 'charset gene_'+str(i+1)+' = DNA, '+ str(seq_length*i+1) + '-' + str(seq_length*i+seq_length) +';\n\t'
partition_info += 'charpartition mine = '
for i in range(len(gene_trees)-1):
    partition_info += 'JC:gene_'+str(i+1)+',\n\t'
partition_info += 'JC:gene_'+str(len(gene_trees))+';\n'
partition_info += 'end;'
print(partition_info)

# Simulate an alignment from the rescaled simulated tree using AliSim
aln, log = simulate_alignment(
        trees = gene_trees,
        partition_info=partition_info,
        partition_type="unlinked",
        rand_seed = seed,
        model='JC',
        num_threads = num_threads,
        length = seq_length,
        population_size = population_size)

# Print the alignment
print(aln)

# Print the console logs
print(log)

```

### Description of Parameters for Tree Simulation
- `samples`: int | dict. The number of sampled individuals in a population.
- `recombination_rate`: float | None. Uniform and non-uniform rates of recombination along the genome.
- `population_size`: int | None. The size of the single constant sized population.
- `ploidy`: int | None. The number of nodes per sample individual, as well as the time scale for continuous time coalescent models. Default: 2, assuming diploid populations.
- `sequence_length`: int | None. The total length of the sequence.
- `model`: str | None. The model under which the ancestral history of the sample is generated. Default: standard coalescent. For other models, please refer to [Discrete Time Wright-Fisher](https://tskit.dev/msprime/docs/stable/api.html#msprime.DiscreteTimeWrightFisher).
- `num_replicates`: int | None. Number of independent simulation replicates (e.g., gene trees).

For further information about the usage of these parameters, see [msprime documentation](https://tskit.dev/msprime/docs/stable/ancestry.html).

### Description of Parameters for Alignment Simulation

- `trees`: list[cogent3.PhyloNode]. A list of cogent3 trees.
- `model`: Model | str. The substitution model. AliSim supports [all models](https://iqtree.github.io/doc/Substitution-Models) available in IQ-TREE. Please refer to [AliSim's User Manual](https://iqtree.github.io/doc/AliSim#specifying-model-parameters) for specifying other substitution models and their parameters.
- `rand_seed`: int. The random seed number.
- `length`: int. The length of sequences. Default: 1000 sites.
- `insertion_rate`: float. The insertion rate. Default: 0.
- `deletion_rate`: float. The deletion rate. Default: 0.
- `insertion_size_distribution`: str. The insertion size distribution. Default: "POW{1.7/100}" (a power-law (Zipfian) distribution with parameter a of 1.7 and maximum indel size of 100).
- `deletion_size_distribution`: str. The deletion size distribution. Default: "POW{1.7/100}".
- `root_seq`: str. The root sequence. Default: "".
- `partition_info`: partition_info: str. The content of the partition file. Default: "".
- `partition_type`: str | None. The type of partitions,  must be ‘equal’, ‘proportion’, ‘unlinked’, or None. Default: None.
- `num_threads`: int. The number of threads. Default: 1.


