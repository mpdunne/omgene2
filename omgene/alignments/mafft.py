import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO
from tempfile import TemporaryDirectory
from typing import List


class Mafft:
    """
    Wrapper class for performing alignments using MAFFT

    """

    def __init__(self, mafft_path: str = '/usr/bin/mafft'):
        """
        Constructor for the MAFFT class.

        :param mafft_path: The path to the MAFFT executable.
        """
        self.mafft_path = mafft_path

    def _run_mafft(self, args: str) -> str:
        """
        Run the MAFFT command.

        :param args: The command line args to pass to MAFFT.
        :return: The result from running MAFFT.
        """
        p = subprocess.Popen(f'{self.mafft_path} --localpair --maxiterate 1000 {args}',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)

        out, err = p.communicate()
        return out.decode()

    def align(self, seqs: List[SeqRecord]) -> List[SeqRecord]:
        """
        Align a list of sequences.

        :param seqs: The sequences.
        :return: The aligned sequences.
        """
        with TemporaryDirectory() as td:
            aa_fasta_file = f'{td}/aa.fa'
            with open(aa_fasta_file, 'w') as f:
                for s in seqs:
                    f.write(f'>{s.id}\n{s.seq}\n')

            result = self._run_mafft(aa_fasta_file)
            return [*SeqIO.parse(StringIO(result), 'fasta')]

    def add(self, aln: List[SeqRecord], seqs: List[SeqRecord]) -> List[SeqRecord]:
        """
        Align sequences to an existing alignment.

        :param aln: The reference alignment.
        :param seqs: The sequences to align to the reference alignment.
        :return: All the aligned sequences.
        """
        with TemporaryDirectory() as td:
            aln_fasta_file = f'{td}/aln.fa'
            with open(aln_fasta_file, 'w') as f:
                for s in aln:
                    f.write(f'>{s.id}\n{s.seq}\n')

            seqs_fasta_file = f'{td}/seqs.fa'
            with open(seqs_fasta_file, 'w') as f:
                for s in seqs:
                    f.write(f'>{s.id}\n{s.seq}\n')

            result = self._run_mafft(f'--add {seqs_fasta_file} {aln_fasta_file}')
            return SeqIO.parse(StringIO(result), 'fasta')
