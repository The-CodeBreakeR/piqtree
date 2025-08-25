# Execute AliSim simulation and output the alignment

One can use `simulate_alignment()` to simulate a multiple sequence alignment, given a list of trees, a model name and a random seed number.

## Usage

### Basic Usage

This example simulates an alignment from a tree under the Juke-Cantor (JC) model and a random seed number of 1. 

```python
from piqtree import simulate_alignment
import cogent3

# Create a cogent3 tree from a Newick string
tree = cogent3.make_tree("((A:0.1,B:0.2):0.1,(C:0.3,D:0.4):0.2,E:0.5);")
trees = [tree]

# Simulate an alignment from the tree(s) under the JC model and a random seed number = 1
aln, log = simulate_alignment(trees, "JC", 1)

# Prints the alignment
print(aln)

# Prints the console logs
print(log)
```

Apart from the simple JC model with no parameters, AliSim also supports all other [more complex models](https://iqtree.github.io/doc/Substitution-Models) available in IQ-TREE. Please refer to [AliSim's User Manual](https://iqtree.github.io/doc/AliSim#specifying-model-parameters) for specifying other substitution models and their parameters.

### Other (Optional) Input Parameters

`simulate_alignment()` also allows specifying the following (optional) parameters.

- `length`: int | list[int]. The length of sequences. Default: 1000 sites. For partitions, if a list of lengths is provided, each length will be applied to the corresponding partition. If a single length is provided, all partitions will have the same length. 
- `insertion_rate`: float. The insertion rate. Default: 0.
- `deletion_rate`: float. The deletion rate. Default: 0.
- `insertion_size_distribution`: str. The insertion size distribution. Default: "POW{1.7/100}" (a power-law (Zipfian) distribution with parameter a of 1.7 and maximum indel size of 100).
- `deletion_size_distribution`: str. The deletion size distribution. Default: "POW{1.7/100}".
- `root_seq`: str. The root sequence. Default: "".
- `partition_type`: str | None. Default: None. If specified, partition_type must be ‘equal’, ‘proportion’, ‘unlinked’, representing Edge-equal, Edge-proportional, and Topology-unlinked [partition models](https://iqtree.github.io/doc/AliSim#partition-models), respectively.  
- `num_threads`: int. The number of threads. Default: 1.
