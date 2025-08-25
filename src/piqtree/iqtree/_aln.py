"""Python wrapper for AliSim Simulation in the IQ-TREE library."""

import tempfile
from pathlib import Path
from typing import Literal

import cogent3
import cogent3.app.typing as c3_types
import yaml
from _piqtree import iq_simulate_alignment
from piqtree.exceptions import ParseIqTreeError

from piqtree.distribution import IndelDistribution
from piqtree.iqtree._decorator import iqtree_func
from piqtree.model import Model

iq_simulate_alignment = iqtree_func(iq_simulate_alignment, hide_files=True)


def simulate_alignment(
        trees: list[cogent3.PhyloNode],
        model: Model | str | list[Model] | list[str],
        rand_seed: int,
        length: int | list[int] = 1000,
        insertion_rate: float = 0.0,
        deletion_rate: float = 0.0,
        insertion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
        deletion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
        root_seq: str | None = None,
        partition_type: Literal["equal", "proportion", "unlinked"] | None = None,
        num_threads: int | None = None,
        population_size: int | None = None,
) -> tuple[c3_types.AlignedSeqsType, str]:
    """Executes AliSim Simulation through IQ-TREE.

    Parameters
    ----------
    trees: list[cogent3.PhyloNode]
        A list of trees.
    subst_model: Model | str | list[Model] | list[str]
        A substitution model. If a list of substitution models is provided, they will serve as partition models.
    rand_seed : int
        The random seed number.
    length: int | list[int], optional
        The length of sequences (by default 1000). If a list of lengths is provided, they will serve as partition lengths.
    insertion_rate: float | None, optional
        The insertion rate (by default 0.0).
    deletion_rate: float | None, optional
        The deletion rate (by default 0.0).
    insertion_size_distribution: str | None, optional
        The insertion size distribution (by default the Zipfian
        distribution with a=1.7 and maximum size 100).
    deletion_size_distribution: str | None, optional
        The deletion size distribution (by default the Zipfian
        distribution with a=1.7 and maximum size 100).
    root_seq: str | None, optional
        The root sequence (by default None).
    partition_type: Literal["equal", "proportion", "unlinked"] | None, optional
        If provided, partition type must be "equal", "proportion", or
        "unlinked" (by default None).
    num_threads: int | None, optional
        Number of threads for IQ-TREE to use, by default None (single-threaded).
    population_size: int | None, optional
        The population size (by default None).

    Returns
    -------
    c3_types.AlignedSeqsType
        The simulated alignment.
    str
        The console log
    """

    if root_seq is None:
        root_seq = ""

    if partition_type is None:
        partition_type = ""

    if num_threads is None:
        num_threads = 1

    if population_size is None:
        population_size = -1

    # convert model and length into lists
    if not isinstance(model, list):
        model = [model]
    if not isinstance(length, list):
        length = [length]

    # check if using partition model, extract the number of partitions
    num_partitions = max(len(trees), len(model), len(length))
    use_partitions = (partition_type is not None and len(partition_type) > 0) or num_partitions > 1

    # if using partitions, partition_type must be specified
    if use_partitions and (partition_type is None or len(partition_type) == 0):
        raise ParseIqTreeError(
            "To use partition models, you must specify partition_type as either `equal` or `proportion`, or `unlinked`.")

    # we don't support `proportion` partition yet, a list of scaling_factors is required
    if partition_type == "proportion":
        raise ParseIqTreeError(
            "Sorry! the current API does not support edge-proportional partitions.")

    # if using partitions, num_partitions must be > 1
    if use_partitions and num_partitions <= 1:
        raise ParseIqTreeError(
            "To use partition models, the number of partitions (i.e., max(len(trees), len(model), len(length)) = " + str(
                num_partitions) + ") must be greater than 1.")

    # the number of models must be either 1 or match the number of partitions
    if use_partitions and len(model) != num_partitions:
        # if one model provided, that model will be used for all partitions
        if len(model) == 1:
            model = model * num_partitions
        # otherwise, throw an error
        else:
            raise ParseIqTreeError(
                "The number of models (" + str(
                    len(model)) + ") must be either 1 or match the number of partitions (i.e., max(len(trees), len(model), len(length)) = " + str(
                    num_partitions) + ").")

    # the number of lengths must be either 1 or match the number of partitions
    if use_partitions and len(length) != num_partitions:
        # if one length provided, that length will be used for all partitions
        if len(length) == 1:
            length = length * num_partitions
        # otherwise, throw an error
        else:
            raise ParseIqTreeError(
                "The number of lengths (" + str(
                    len(length)) + ") must be either 1 or match the number of partitions (i.e., max(len(trees), len(model), len(length)) = " + str(
                    num_partitions) + ").")

    # convert model to string if needed
    for i in range(len(model)):
        model[i] = str(model[i]) if isinstance(model[i], Model) else model[i]

    # Convert the trees to Newick strings
    newick_trees = [str(tree) for tree in trees]

    # init parameters for AliSim API
    common_model = model[0]
    common_length = length[0]

    # generate a partition file if needed
    partition_info = ""
    if use_partitions:
        partition_info = '#nexus\n begin sets;\n\t'
        start_site = 0
        for i in range(num_partitions):
            partition_info += 'charset gene_' + str(i + 1) + ' = ' + str(start_site + 1) + '-' + str(
                start_site + length[i]) + ';\n\t'
            start_site += length[i]
        partition_info += 'charpartition mine = '
        for i in range(num_partitions - 1):
            partition_info += model[i] + ':gene_' + str(i + 1) + ',\n\t'
        partition_info += model[-1] + ':gene_' + str(num_partitions) + ';\n'
        partition_info += 'end;'

    # Call the IQ-TREE function
    yaml_result = yaml.safe_load(
        iq_simulate_alignment(
            newick_trees,
            common_model,
            rand_seed,
            partition_info,
            partition_type,
            common_length,
            insertion_rate,
            deletion_rate,
            root_seq,
            num_threads,
            str(insertion_size_distribution),
            str(deletion_size_distribution),
            population_size,
        ),
    )

    # Extract the simulated alignment and the content of the log file from the YAML result
    # cogent3.make_aligned_seqs dictionary
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
        tmp.write(yaml_result["alignment"])
        tmp.flush()
        aln = cogent3.load_aligned_seqs(
            tmp.name,
            format="fasta",
            new_type=True,
        )  # moltype = 'dna' or 'protein')
    Path(tmp.name).unlink()
    console_log = yaml_result["log"]

    return aln, console_log
