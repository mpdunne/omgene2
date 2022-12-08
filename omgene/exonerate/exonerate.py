import subprocess
import tempfile
import pandas as pd

from typing import List, Union

exonerate_path = '~/software/exonerate-2.2.0-x86_64/bin/exonerate'


class Exonerate:
    """
    Wrapper class for the Exonerate package, to simply run exonerate and parse its output.

    """

    def __init__(self, exonerate_path: str = ''):
        """
        Constructor for the Exonerate class.

        :param exonerate_path: The path to the exonerate executable
        """
        self.exonerate_path = exonerate_path
        self._check_exonerate_path()

    def _check_exonerate_path(self):
        """
        Runs the help method of exonerate to check that it can be run.

        :return: Nothing if everything's fine. Raise error if not.
        """
        p = subprocess.Popen(f'{self.exonerate_path} -h',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)

        out, err = p.communicate()
        assert not err
        assert 'exonerate from exonerate' in out.decode()

    def _run_exonerate(self, aa_fasta_file: str, genome_fasta_file: str) -> str:
        """
        Run a basic exonerate command. Searches the provided genome_fasta_file using the
        sequences in aa_fasta_file.

        :param aa_fasta_file: The file containing the AA sequences.
        :param genome_fasta_file: The file containing the genome sequences to search.
        :return: The output from exonerate, in bytes format.
        """
        p = subprocess.Popen(f'{self.exonerate_path} --showalignment no '
                             f'--showtargetgff --frameshift -10000 --model protein2genome '
                             f'{aa_fasta_file} {genome_fasta_file}',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)

        out, err = p.communicate()
        return out.decode()

    def _process_output(self, output: str) -> List[pd.DataFrame]:
        """
        Process the exonerate output straight from the command line.

        :param output: The output from the exonerate program in bytes.
        :return: A list of LOLs containing GFF information.
        """
        all_gffs = []
        current_gff = []
        store_lines = False

        for line in output.split('\n'):
            if line == '# --- START OF GFF DUMP ---':
                store_lines = True
            elif line == '# --- END OF GFF DUMP ---':
                store_lines = False

                current_gff_df = pd.DataFrame(current_gff)
                current_gff_df.columns = ['chr', 'source', 'type', 'start', 'end', '.', 'strand', 'frame', 'desc']

                all_gffs.append(current_gff_df)
                current_gff = []

            if not store_lines or line[0] == '#':
                continue
            else:
                current_gff.append(line.split('\t'))

        return all_gffs

    def run(self, aa_sequences: Union[List[str], str], genome_sequence: str) -> List[pd.DataFrame]:
        """
        Do a basic run of exonerate and process the output.

        :param aa_sequences: The AA sequences (or single sequence) use.
        :param genome_sequence: The DNA sequence to search.
        :return: A list of GFFs in list-of-lists format.
        """
        if type(aa_sequences) is str:
            aa_sequences = [aa_sequences]

        with tempfile.TemporaryDirectory() as td:
            aa_fasta_file = f'{td}/aa.fa'
            genome_fasta_file = f'{td}/genome.fa'

            with open(aa_fasta_file, 'w') as f:
                for i, seq in enumerate(aa_sequences):
                    f.write(f'>seq_{i}\n{seq}\n')

            with open(genome_fasta_file, 'w') as f:
                f.write(f'>genome\n{genome_sequence}\n')

            output = self._run_exonerate(aa_fasta_file, genome_fasta_file)
            gffs = self._process_output(output)

            return gffs
