from multiprocessing import Pool, cpu_count
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio.Seq import Seq

def build_k_mers(reads: list[str], k: int) -> list[list[str]]:
  """
  Builds and returns all k-mers from the given reads.

  :param reads: Reads.
  :type reads: list[str]
  :param k: Length of the k-mers. All reads shorter than k are discarted.
  :type k: int

  :return: List of lists, k-mers for each read.
  :rtype: list[list[str]]
  """
  return [[read[i : i + k] for i in range(len(read) - k + 1)] for read in reads if len(read) >= k]

def build_k_mers_helper(args: tuple[str, int]) -> list[str]:
  """
  Helper method for build_k_mers_parallel.
  Builds the list of k-mers of a single sequence.

  :param args: Tuple with the sequence and the length of the k-mers.
  :type args: tuple[str, int]

  :return: List of all k-mers for the given sequence.
  :rtype: list[str]
  """
  seq, k = args
  k_mers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
  return k_mers

def build_k_mers_parallel(reads: list[str], k: int) -> list[list[str]]:
  """
  Builds and returns all k-mers from the given reads.
  Runs in parallel.

  :param reads: Reads.
  :type reads: list[str]
  :param k: Length of the k-mers. All reads shorter than k are discarted.
  :type k: int

  :return: List of lists, k-mers for each read.
  :rtype: list[list[str]]
  """
  with Pool(cpu_count()) as pool:
    k_mers = pool.map(build_k_mers_helper, [(read, k) for read in reads if len(read) >= k])
  return k_mers

def build_k_mers_helper_batch(args: tuple[list[str], int]) -> list[list[str]]:
  """
  Helper method for build_k_mers_parallel_batch.
  Builds the list of k-mers of a single batch of sequences.

  :param args: Tuple with the batch of sequences and the length of the k-mers.
  :type args: tuple[list[str], int]

  :return: List of all k-mers for the given batch.
  :rtype: list[str]
  """
  batch, k = args
  return [[seq[i : i + k] for i in range(len(seq) - k + 1)] for seq in batch]

def build_k_mers_parallel_batch(reads: list[str], k: int) -> list[list[str]]:
  """
  Builds and returns all k-mers from the given reads.
  Runs in parallel with batches.

  :param reads: Reads.
  :type reads: list[str]
  :param k: Length of the k-mers. All reads shorter than k are discarted.
  :type k: int

  :return: List of lists, k-mers for each read.
  :rtype: list[list[str]]
  """
  seqs = [read for read in reads if len(read) >= k]
  step = int(len(seqs) / cpu_count()) or 1
  with Pool(cpu_count()) as pool:
    k_mers_batches = pool.map(build_k_mers_helper_batch, [(seqs[i : i + step], k) for i in range(0, len(seqs), step)])
  return [k_mers for batch in k_mers_batches for k_mers in batch]
