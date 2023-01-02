import re
import numpy as np

from omgene.alignments.mafft import Mafft
from omgene.alignments.score import MSAScorer
from omgene.alignments.match import get_nonmatching_regions
from omgene.exonerate.exonerate import Exonerate
from omgene.utils.maps import *
from omgene.utils.gff import GeneContext
from omgene.utils.features import find_start, find_stop

from copy import deepcopy
from typing import Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from typing import List, Tuple


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

    def run(self, return_scores: bool = False) -> Dict[str, GeneContext]:
        """
        Run the optimisation.

        :param return_scores: Whether or not to return the old and new MSA scores.
        :return: Optimised genes, keyed by transcript ID.
        """
        print('Running search...')
        options = self._exonerate_all_v_all()

        print('Adjusting stop codons')
        options = {k: [find_stop(o) for o in opts] for k, opts in options.items()}
        options = {k: [o for o in opts if o is not None] for k, opts in options.items()}

        print('Adding start codons')
        valid_genes = {k: [] for k in options.keys()}
        for k, opts in options.items():
            for opt in opts:
                if opt.seq[:3] == 'ATG':
                    valid_genes[k].append(opt)
                else:
                    left_option = find_start(opt, direction='left')
                    if left_option is not None:
                        valid_genes[k].append(left_option)

                    right_option = find_start(opt, direction='right')
                    if right_option is not None:
                        valid_genes[k].append(right_option)

        print('Preparing sequences for optimisation...')
        original_tids, all_gcs, all_aligned_seqs = self._prepare_candidates(valid_genes)

        print('Choosing best sequences...')
        result = self._optimise_gene_models(original_tids, all_gcs, all_aligned_seqs, return_scores=return_scores)
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

    def _prepare_candidates(self, candidates: Dict[str, List[GeneContext]]) -> \
        Tuple[Dict[str, str], Dict[str, Dict[str, GeneContext]], Dict[str, Seq]]:
        """
        Organise and align the candidates in order for them to be processed.

        Returns:
            - A dict specifying the original transcript id for each gene id.
            - A dict, keyed by transcript ID, containing the GeneContext for each candidate transcript
            - A dict, keyed by gene ID then transcript ID, containing the aligned sequence of each candidate transcript.

        :param candidates: The candidates as returned by the cross-species search.
        :return: A tuple of (original_tids, gene_contexts, aligned_seqs)
        """
        # Organise the new candidates and the original sequences.
        all_gcs = {}
        original_tids = {}

        for gid, original_gc in {**self.reference_genes, **self.target_genes}.items():
            all_gcs[gid] = {gid: original_gc}
            original_tids[gid] = gid

        for gid, candidate_gcs in candidates.items():
            candidate_gcs_by_tid = {f'{gid}.alt-{i}': gc for i, gc in enumerate(candidate_gcs)}
            all_gcs[gid].update(candidate_gcs_by_tid)

        # Now we want to align the originals and then align all the new sequences to them.
        # First get all original seqs and all options in list format.
        seq_list_orig = []
        seq_list_candidates = []

        for gid, original_tid in original_tids.items():

            # Original sequence
            orig_gc = all_gcs[gid][original_tid]
            seq_list_orig.append(SeqRecord(Seq(orig_gc.seq.translate()), id=gid))

            # The candidates are then anything that's not the original
            for tid, candidate_gc in all_gcs[gid].items():
                if tid != original_tid:
                    seq_list_candidates.append(SeqRecord(Seq(all_gcs[gid][tid].seq.translate()), id=tid))

        # Run alignments on the original sequences as well as all the options....
        mafft = Mafft()

        print('Aligning original sequences...')
        aln_orig = [*mafft.align(seq_list_orig)]

        print('Aligning new candidates...')
        aln_new = [*mafft.add(aln_orig, seq_list_candidates)]
        aln_new_seqs = {s.id: str(s.seq) for s in aln_new}

        return original_tids, all_gcs, aln_new_seqs

    def _optimise_gene_models(self,
                              original_tids: Dict[str, str],
                              all_gcs: Dict[str, Dict[str, GeneContext]],
                              all_aligned_seqs: Dict[str, Seq],
                              return_scores=True):
        """
        Choose the best ones.

        :param original_tids: A dict matching each gene ID to its original transcipt ID.
        :param all_gcs: The GeneContexts for all the candidates.
        :param all_aligned_seqs: A dict of dicts containing all the aligned candidate sequences.
        :param return_scores: Whether to return the scores.
        :return: The best option for each transcript ID.
        """

        print('Scoring and choosing...')
        m = MSAScorer()

        # Our starting point
        current_msa = {gid: all_aligned_seqs[tid] for gid, tid in original_tids.items()}
        current_gcs = {gid: all_gcs[gid][tid] for gid, tid in original_tids.items()}
        original_score = m.alignment_score([*current_msa.values()])

        # Now go through each sequence and greedily look for improvements.
        repeat = True
        while repeat:
            repeat = False
            for gid, alternatives in all_gcs.items():
                for tid_alt, gc_alt in alternatives.items():
                    if current_gcs[gid].seq.translate() == gc_alt.seq.translate():
                        continue

                    # Unpack the aligned sequences for original and alt.
                    aligned_seq_current = current_msa[gid]
                    aligned_seq_alt = all_aligned_seqs[tid_alt]

                    nonmatching_regions = get_nonmatching_regions(aligned_seq_current, str(aligned_seq_alt))
                    if not nonmatching_regions:
                        continue

                    # Get score differentials for nonmatching regions.
                    nonmatching_region_score_differentials = []

                    alt_msa = {**current_msa, **{gid: aligned_seq_alt}}
                    current_msa_array = np.array([list(seq) for seq in current_msa.values()])
                    alt_msa_array = np.array([list(seq) for seq in alt_msa.values()])

                    for (l, r) in nonmatching_regions:
                        score_current = np.sum(m.column_scores(current_msa_array[:, l: r]))
                        score_alt = np.sum(m.column_scores(alt_msa_array[:, l: r]))
                        nonmatching_region_score_differentials.append(score_alt - score_current)

                    best_ix = np.argmax(nonmatching_region_score_differentials)
                    best_region_candidate = None

                    if nonmatching_region_score_differentials[best_ix] > 0:
                        best_region_candidate = nonmatching_regions[best_ix]

                    if best_region_candidate is None:
                        continue

                    # Get GCs and exon maps.
                    gc_orig = current_gcs[gid]
                    gc_alt = all_gcs[gid][tid_alt]

                    map_orig = get_exon_map(gc_orig)
                    map_alt = get_exon_map(gc_alt)

                    # Get AA maps
                    aa_map_orig = exon_map_to_aa_map(map_orig)
                    aa_map_alt = exon_map_to_aa_map(map_alt)

                    gene_map_region_alt = aligned_seq_region_to_exon_map_region(aligned_seq_alt, best_region_candidate,
                                                                                aa_map_alt)

                    # Attempt a transplant
                    try:
                        exon_map_transplanted = exon_map_transplant(map_orig, map_alt, gene_map_region_alt)
                    except TransplantError:
                        continue

                    # Get the exons back
                    exons_transplanted = exons_from_exon_map(exon_map_transplanted)
                    gc_transplanted = deepcopy(gc_orig)
                    gc_transplanted.meta_exons = exons_transplanted

                    # Update current
                    current_gcs[gid] = gc_transplanted
                    l, r = best_region_candidate
                    current_aligned_seq = aligned_seq_current[:l] + aligned_seq_alt[l:r] + aligned_seq_current[r:]
                    current_msa[gid] = current_aligned_seq

                    # If we've got here, do it all again.
                    repeat = True

        current_score = m.alignment_score([*current_msa.values()])
        if current_score > original_score:
            print(f'Improved MSA score from {original_score:4f} to {current_score:4f}.')
        else:
            print('No improvement.')

        if return_scores:
            return current_gcs, original_score, current_score
        else:
            return current_gcs
