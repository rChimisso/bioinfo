from multiprocessing import Pool, cpu_count
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio.Seq import Seq

def build_k_mers(reads: FastqPhredIterator, k: int) -> list[Seq]:
  """
  Builds and returns all k-mers from the given list of reads.
  Instead of a list of k-mers for each read, the returned list is flattened.
  Runs in series.

  :param reads: Iterator of reads.
  :type reads: FastqPhredIterator
  :param k: Length of the k-mers. All reads shorter than k are discarted.
  :type k: int

  :return: Flattened list of all k-mers for each read.
  :rtype: list[Seq]
  """
  seqs: list[Seq] = [read.seq for read in reads if len(read) >= k]
  print(f'Number of reads that can be used to build the k-mers: {len(seqs)}')
  k_mers = [[seq[i : i + k] for i in range(len(seq) - k + 1)] for seq in seqs]
  print('Terminated building all k-mers, proceeding with flattening.')
  return [k_mer for seq_k_mers in k_mers for k_mer in seq_k_mers]

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

def build_k_mers_parallel(reads: FastqPhredIterator, k: int) -> list[Seq]:
  """
  Builds and returns all k-mers from the given list of reads.
  Instead of a list of k-mers for each read, the returned list is flattened.
  Runs in parallel.

  :param reads: Iterator of reads.
  :type reads: FastqPhredIterator
  :param k: Length of the k-mers. All reads shorter than k are discarted.
  :type k: int

  :return: Flattened list of all k-mers for each read.
  :rtype: list[Seq]
  """
  seqs: list[Seq] = [read.seq for read in reads if len(read) >= k]
  print(f'Number of reads that can be used to build the k-mers: {len(seqs)}')
  with Pool(cpu_count()) as pool:
    k_mers = pool.map(build_k_mers_helper, [(seq, k) for seq in seqs])
  print('Terminated building all k-mers, proceeding with flattening.')
  return [k_mer for seq_k_mers in k_mers for k_mer in seq_k_mers]
