from Bio import SeqIO
from Bio.Seq import Seq

def export_fastq(input: str, output: str) -> None:
  """
  Exports the reads from the specified FASTQ file into the specified TXT file.

  :param input: Name of the FASTQ file (without extension).
  :type input: str
  :param output: Name of the TXT file (without extension).
  :type output: str
  """
  with open(f'{input}.fastq') as fastq:
    with open(f'{output}.txt', 'w+') as file:
      file.writelines([str(seq) + '\n' for read in SeqIO.parse(fastq, "fastq") if len((seq := read.seq)) > 0])

def extract_seqs(input: str) -> list[Seq]:
  """
  Returns the sequences contained in the specified TXT file, wrapping them as Seq objects.

  :param input: Name of the TXT file (without extension).
  :type input: str
  """
  with open(f'{input}.txt') as file:
    return [Seq(seq) for seq in file.readlines()]
