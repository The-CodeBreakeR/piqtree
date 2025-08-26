import pytest
import cogent3
from cogent3.core.alignment import Alignment

from piqtree import simulate_alignment
from piqtree.exceptions import ParseIqTreeError

trees = [
        cogent3.make_tree("((A:0.1,B:0.2):0.1,(C:0.3,D:0.4):0.2,E:0.5);"),
        cogent3.make_tree("(((F:0.2,E:0.4):0.3,(C:0.1,B:0.5):0.2):0.4,(A:0.3,D:0.2):0.5);"),
    ]

@pytest.mark.parametrize("tree", trees)
@pytest.mark.parametrize("model", ["JC", "HKY{2}+F{0.2/0.3/0.1/0.4}"])
@pytest.mark.parametrize("rand_seed", [1, 3])
@pytest.mark.parametrize("length", [100, 2000])
@pytest.mark.parametrize("insertion_rate", [0, 0.02])
@pytest.mark.parametrize("deletion_rate", [0, 0.04])
@pytest.mark.parametrize("insertion_size_distribution", ["", "NB{5/20}"])
@pytest.mark.parametrize("deletion_size_distribution", ["POW{1.7/100}"])
@pytest.mark.parametrize("num_threads", [None, 4])
def test_simulate_alignment_basic(
   tree: cogent3.PhyloNode,
   model: str,
   rand_seed: int,
   length: int | None,
   insertion_rate: float | None,
   deletion_rate: float | None,
   insertion_size_distribution: str,
   deletion_size_distribution: str,
   num_threads: int | None,
) -> None:
   res = simulate_alignment(
       tree = tree,
       model = model,
       rand_seed = rand_seed,
       length = length,
       insertion_rate = insertion_rate,
       deletion_rate = deletion_rate,
       insertion_size_distribution = insertion_size_distribution,
       deletion_size_distribution = deletion_size_distribution,
       num_threads = num_threads,
   )


   # Checks if simulate_alignment is returning a tuple with two values
   assert isinstance(res, tuple) and len(res) == 2


   aln, log = res


   # Checks if the number of sequences in the output alignment is equal to the number of taxa in the input tree
   assert aln.num_seqs == len(tree.tips())


   # If no insertion or deletion rates are given, checks if every sequence length in the output alignment is equal to the input sequence length
   if insertion_rate is None and deletion_rate is None:
       unique_seq_lengths = list(set(aln.get_lengths().to_dict().values()))
       actual_seq_length = seq_length if seq_length is not None else 1000
       assert len(unique_seq_lengths) == 1 and unique_seq_lengths[0] == actual_seq_length



# ---------------------------
# (MORE ADVANCED) SUCCESSFUL CASES
# ---------------------------


def test_single_model_multiple_trees():
   aln, log = simulate_alignment(tree=trees, model="HKY{2}+F{0.2/0.3/0.1/0.4}", rand_seed=123, length=[10, 20], partition_type="unlinked")
   assert "ACGT" in str(aln.to_dict())
   assert isinstance(log, str)




def test_multiple_models():
   models = ["JC", "GTR{2/3/4/5/6}"]
   lengths = [50, 100]
   aln, log = simulate_alignment(tree=trees, model=models, rand_seed=11, length=lengths, partition_type="unlinked")
   assert aln is not None




def test_with_population_size():
   aln, log = simulate_alignment(tree=trees[0], model="JC", rand_seed=99, population_size=100)
   assert aln is not None




def test_with_root_seq():
   aln, log = simulate_alignment(tree=trees[0], model="GTR{2/3/4/5/6}+F{0.2/0.3/0.1/0.4}", rand_seed=1, root_seq="ACGT")
   assert aln is not None


# ---------------------------
# ERROR CASES
# ---------------------------


def test_invalid_population_size():
   with pytest.raises(ParseIqTreeError, match="population_size must be positive"):
       simulate_alignment(tree=trees[0], model="HKY{2}", rand_seed=1, population_size=0)




def test_partition_type_required():
   with pytest.raises(ParseIqTreeError, match="must specify partition_type"):
       simulate_alignment(tree=trees, model="HKY{3}+F{0.2/0.3/0.1/0.4}", rand_seed=2, length=[10, 20])




def test_proportion_partition_not_supported():
   with pytest.raises(ParseIqTreeError, match="does not support edge-proportional"):
       simulate_alignment(tree=trees, model=["JC", "HKY{2}"], rand_seed=3,
                          length=[10, 20], partition_type="proportion")




def test_equal_partition_required_one_tree():
   with pytest.raises(ParseIqTreeError, match="only one tree is required"):
       simulate_alignment(tree=trees, model=["HKY{3}+F{0.2/0.3/0.1/0.4}", "JC"], rand_seed=3,
                          length=[10, 20], partition_type="equal")




def test_mismatched_tree_count():
   with pytest.raises(ParseIqTreeError, match="number of trees"):
       simulate_alignment(tree=trees, model=["JC"], rand_seed=4,
                          length=[10, 20, 30], partition_type="unlinked")




def test_mismatched_model_count():
   with pytest.raises(ParseIqTreeError, match="number of models"):
       simulate_alignment(
           tree=trees.extend(trees),
           model=["HKY{3}+F{0.2/0.3/0.1/0.4}", "GTR{2/3/4/5/6}+F{0.2/0.3/0.1/0.4}", "JC"],
           rand_seed=5,
           length=[10, 20, 30, 40],
           partition_type="unlinked"
       )




def test_mismatched_length_count():
   with pytest.raises(ParseIqTreeError, match="number of lengths"):
       simulate_alignment(
           tree=trees.extend(trees[0]),
           model=["HKY{3}+F{0.2/0.3/0.1/0.4}", "GTR{2/3/4/5/6}+F{0.2/0.3/0.1/0.4}", "JC"],
           rand_seed=6,
           length=[10, 20],
           partition_type="unlinked"
       )

# WRONG DATA TYPES

def test_tree_wrong_type():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree="not_a_tree", model="JC", rand_seed=1)
    print("Caught error message:", excinfo.value)

def test_model_wrong_type():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree=trees[0], model=12345, rand_seed=1)
    print("Caught error message:", excinfo.value)

def test_length_wrong_type():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree=trees[0], model="JC", rand_seed=1, length="1000")
    print("Caught error message:", excinfo.value)

def test_partition_type_wrong_type():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree=trees, model="JC", rand_seed=1, partition_type=123)
    print("Caught error message:", excinfo.value)

# OUT OF RANGE VALUES

def test_negative_length():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree=trees[0], model="JC", rand_seed=1, length=-100)
    print("Caught error message:", excinfo.value)

def test_negative_insertion_rate():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree=trees[0], model="JC", rand_seed=1, insertion_rate=-0.01)
    print("Caught error message:", excinfo.value)

def test_negative_deletion_rate():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree=trees[0], model="JC", rand_seed=1, deletion_rate=-0.01)
    print("Caught error message:", excinfo.value)

def test_population_size_negative():
    with pytest.raises(Exception) as excinfo:
        simulate_alignment(tree=trees[0], model="JC", rand_seed=1, population_size=-10)
    print("Caught error message:", excinfo.value)

