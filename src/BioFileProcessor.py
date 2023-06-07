from typing import Generator, Iterator, Union
from Bio import SeqIO
from pathlib import Path
import log as log

class BioFileProcessor:
  """
  Processor for biological reads data files in FASTA or FASTQ format.

  :ivar filepath: Path to the data file.
  :vartype filepath: str
  :ivar format: Format of the data file, either 'fasta' or 'fastq'.
  :vartype format: str
  """

  def __init__(self, filepath: str):
    """
    Constructs all the necessary attributes for the BioFileProcessor object.

    :param filepath: The path to the data file.
    :type filepath: str
    :raises ValueError: If the file format is neither 'fasta' nor 'fastq'.
    """
    self.filepath: str = filepath
    self.format: Union[str, None] = None
    if Path(filepath).suffix == '.fasta':
      self.format = 'fasta'
    elif Path(filepath).suffix == '.fastq':
      self.format = 'fastq'
    else:
      raise ValueError(f'File extension not recognized for {filepath}. Use .fasta or .fastq files.')
    log.logger().info(f'BioFileProcessor for the file [{filepath}] created successfully') 

  def is_acgt_read(self, read: str) -> bool:
    """
    Checks if a read is composed only of the characters 'a', 'c', 'g', 't'.

    :param read: DNA read to check.
    :type read: str
    :return: True if the read is valid, False otherwise.
    :rtype: bool
    """
    return all(base in 'acgt' for base in read.lower())

  def iter_filtered_reads(self, k: int) -> Iterator[str]:
    """
    Generates valid reads from the data file one at a time.

    :param k: K-mer length.
    :type k: int
    :return: An Iterator object that yields valid reads one at a time.
    :rtype: Iterator[str]
    """
    for read in SeqIO.parse(self.filepath, self.format):
      if len((seq := str(read.seq))) >= k and self.is_acgt_read(read.seq):
        yield seq
