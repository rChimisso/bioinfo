from multiprocessing import Pool, cpu_count
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio.Seq import Seq

def build_k_mers(reads: FastqPhredIterator, k: int) -> list[list[Seq]]:
  """
  Builds and returns all k-mers from the given reads.

  :param reads: Iterator of reads.
  :type reads: FastqPhredIterator
  :param k: Length of the k-mers. All reads shorter than k are discarted.
  :type k: int

  :return: List of lists, k-mers for each read.
  :rtype: list[Seq]
  """
  return [[seq[i : i + k] for i in range(len(seq) - k + 1)] for read in reads if len((seq := read.seq)) >= k]

def build_k_mers_helper(args: tuple[Seq, int]) -> list[Seq]:
  """
  Helper method for build_k_mers_parallel.
  Builds the list of k-mers of a single sequence.

  :param args: Tuple with the sequence and the length of the k-mers.
  :type args: tuple[Seq, int]

  :return: List of all k-mers for the given sequence.
  :rtype: list[Seq]
  """
  seq, k = args
  k_mers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
  return k_mers

def build_k_mers_parallel(reads: FastqPhredIterator, k: int) -> list[list[Seq]]:
  """
  Builds and returns all k-mers from the given reads.
  Runs in parallel.

  :param reads: Iterator of reads.
  :type reads: FastqPhredIterator
  :param k: Length of the k-mers. All reads shorter than k are discarted.
  :type k: int

  :return: List of lists, k-mers for each read.
  :rtype: list[Seq]
  """
  with Pool(cpu_count()) as pool:
    k_mers = pool.map(build_k_mers_helper, [(seq, k) for read in reads if len((seq := read.seq)) >= k])
  return k_mers
