from omgene.alignments.mafft import Mafft
from omgene.alignments.score import MSAScorer
from omgene.exonerate.exonerate import Exonerate
from omgene.utils.gff import GeneContext

from copy import deepcopy
from typing import Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from typing import List


class OMGene:
    """
    Main class for running the OMGene algorithm.

    """

    def __init__(self,
                 exonerate_path='~/software/exonerate-2.2.0-x86_64/bin/exonerate',
                 mafft_path='/usr/bin/mafft'):
        """
        Constructor for the OMGene class.

        """
        self.reference_genes = {}
        self.target_genes = {}
        self.exonerate_path = exonerate_path
        self.mafft_path = mafft_path

    def add_genes(self, genes: Dict[str, GeneContext], reference: bool = False):
        """
        Add genes. Reference genes will be used for search and comparison but will not be optimised.

        :param genes: The genes to add.
        :param reference: Whether we are adding them as reference genes or not.
        :return: None.
        """
        if reference:
            self.reference_genes.update(genes)
        else:
            self.target_genes.update(genes)

    def run(self) -> Dict[str, GeneContext]:
        """
        Run the optimisation.

        :return: Optimised genes, keyed by transcript ID.
        """
        print('Running search...')
        options = self._exonerate_all_v_all()
        valid_genes = {k: [o for o in opts if o.seq[0:3] == 'ATG'] for k, opts in options.items()}

        print('Choosing best sequences...')
        result, original_score, new_score = self._choose_best(valid_genes)
        # print(original_score, new_score)
        return result

    def _exonerate_all_v_all(self) -> Dict[str, List[GeneContext]]:
        """
        Cross-run exonerate on all of the stored genes.

        :return: A dict of options for each transcript ID key.
        """

        e = Exonerate(self.exonerate_path)

        query_genes = {**self.reference_genes, **self.target_genes}
        target_genes = self.target_genes

        options = {tid: [] for tid in target_genes.keys()}

        for target_tid, target_data in target_genes.items():
            print(f'Searching against {target_tid}')

            target_genome_seq = target_data.meta_chr_seq
            query_seqs = [query_data.seq.translate() for query_data in query_genes.values()]

            results = e.run(query_seqs, target_genome_seq)
            for result in results:
                new_target_data = deepcopy(target_data)
                result_cds = result[result['type'] == 'cds']

                new_target_data.meta_exons = [[int(row['start']) - 1, int(row['end'])] for _, row in
                                              result_cds.iterrows()]
                options[target_tid].append(new_target_data)

        return options

    def _choose_best(self, options: Dict[str, List[GeneContext]], return_scores=True):
        """
        Choose the best ones.

        :param options: The options for each gene, keyed by transcript ID.
        :param return_scores: Whether to return the scores.
        :return: The best option for each transcript ID.
        """
        options_named = {s: {f'{s}.alt-{i}': o for i, o in enumerate(os)} for s, os in options.items()}

        query_genes = {**self.reference_genes, **self.target_genes}
        seqs_orig = [SeqRecord(Seq(o.seq.translate()), id=s) for s, o in query_genes.items()]
        seqs_opts = [SeqRecord(Seq(o.seq.translate()), id=s) for _, os in options_named.items() for s, o in os.items()]

        mafft = Mafft()

        print('Aligning original sequences...')
        aln = [*mafft.align(seqs_orig)]

        print('Aligning new candidates...')
        aln_new = [*mafft.add(aln, seqs_opts)]
        aln_new_seqs = {s.id: str(s.seq) for s in aln_new}

        print('Scoring and choosing...')
        m = MSAScorer()
        current = {s: s for s in query_genes.keys()}
        potentials = {s: [f'{s}.alt-{i}' for i, o in enumerate(os)] for s, os in options.items()}
        current_score = original_score = m.alignment_score([aln_new_seqs[k] for k in current.values()])

        best = current
        best_score = current_score
        again = True

        while again:
            again = False
            current = best.copy()
            for tid, new_tids in potentials.items():
                for new_tid in new_tids:
                    current[tid] = new_tid
                    current_score = m.alignment_score([aln_new_seqs[k] for k in current.values()])
                    if current_score > best_score:
                        again = True
                        best = current.copy()
                        best_score = current_score

        result = {}
        for s_orig, s_new in best.items():
            if s_orig == s_new:
                result[s_orig] = query_genes[s_orig]
            else:
                result[s_orig] = options_named[s_orig][s_new]

        if best_score > original_score:
            print(f'Improved MSA score from {original_score:4f} to {best_score:4f}.')
        else:
            print('No improvement.')

        if return_scores:
            return result, original_score, best_score
        else:
            return result
