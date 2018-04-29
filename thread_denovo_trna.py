#!/usr/bin/env python3

import sys, subprocess, mypy
from typing import Dict, List, Tuple
from lib.Sequence import Sequence
from lib.file_io import exe, my_loc, get_seqs, get_seqs_mfa
from lib.str_manip import modomics_from_pdb, U_equivs, A_equivs, C_equivs, G_equivs
from lib.rosetta_calls import remodel_new_sequence
from timeit import timeit, repeat

"""
Homology modeling of novel tRNA sequences expressed in Modomics's extended
single-letter-code nomenclature.

* Distributed with this software *
 - a multiple sequence alignment of every atomic resolution tRNA in the PDB 
   (omitting massive ribosome structures)
 - those structures, trimmed and renormalized to make them robust for Rosetta

* Your responsibility *
 - Provide a tRNA sequence. You can manually align it against our MSA or we can
   do it for you.
 - Select a modeling method. Your options:
   + FARFAR -- relatively fast and relatively low-resolution, should yield 
     ~10 structures per CPU-hour
   + STEPWISE -- relatively slow but very high-resolution, should yield 
     1 structure per ~6 CPU-hours. (Use for structures with a few deletions 
     making key contacts, like if the T-loop residues G18-G19 are altered vs.
     the best available template.)
 - Provide a number of desired structures to produce. This many structures will
   be produced for *every* template that's identified.

 * Optional *
  - If there is a structure you believe ought to be a template that isn't coming
    up (perhaps it's your structure, which you haven't deposited yet!), you can 
    manually add a structure to the template library with its sequence alignment
    (or edit an existing alignment).
  - If no sequence in the PDB is a good template match to your target sequence,
    then you can instead operate in 'density mode', which attempts template-free
    modeling guided instead by weak electron density and base pairing 
    constraints.
"""


###
#   Sequence alignment scoring
###

import numpy as np

def simple_match(the_seq1: Sequence, the_seq2: Sequence, quiet=True) -> int:
    """
    Assumes seq1 and seq2 are aligned already (i.e., a 'static score').
    Rewards similarity in length, letter matches, purine matches, and
    variant matches.

    Ooh, consider not penalizing complementary WC mismatches within helices?
    would need SS. But, now we have that! So we can turn this on whenever we want.

    NOTE: I have tried using numba's JIT as well as making seq1, seq2 into
    numpy objects. The former didn't run and the latter made it 2-5x slower.
    This is still THE bottleneck though.
    """

    seq1: str = the_seq1.sequence
    seq2: str = the_seq2.sequence

    score: int = (len(seq1) - len(seq2))*3
    
    # For numpy chararrays, must iterate via index.

    for c1, c2 in zip(seq1, seq2):
        if c1 == c2: score += 5
        else:
            if c1 == 'A' and c2 == 'G': score += 2
            elif c1 == 'G' and c2 == 'A': score += 2
            elif c1 == 'C' and c2 == 'U': score += 2
            elif c1 == 'U' and c2 == 'C': score += 2

            elif c1 == 'U' and c2 in U_equivs: score += 3
            elif c1 in U_equivs and c2 == 'U': score += 3
            elif c1 == 'A' and c2 in A_equivs: score += 3
            elif c1 in A_equivs and c2 == 'A': score += 3
            elif c1 == 'G' and c2 in G_equivs: score += 3
            elif c1 in G_equivs and c2 == 'G': score += 3
            elif c1 == 'C' and c2 in C_equivs: score += 3
            elif c1 in C_equivs and c2 == 'C': score += 3

            elif (c1 == '-' or c2 == '-'): score -= 10

    if not quiet:
        pass
        #print(score)
        #print(seq1)
        #print(seq2)
    return score


"""
DOMAIN ANALYSIS:

1. Proximal stuff
2. Never rebuild key tert contact (GG t-loop) unless it is directly altered 
"""

def match_seq_to_best_template_seq(tgt_seq: Sequence, templates: List[Tuple[str, Sequence]]) -> List[Tuple[str, Sequence]]:
    print("\nMatching desired target sequence to our template library...")

    # Filter sequences for those of same true length. If 1, done.
    # If zero, quit entirely; no point.
    new_seqs: List[Tuple[str, Sequence]] = [t for t in templates]
    ii: int = 0

    pdb_len: int = len(new_seqs[0][1].sequence.replace('-',''))

    # HOLY FUCK TRY SOME ALIGNMENTS YOU MORON

    seq_matching_to_pdb: Dict[str, int] = {}
    for ii, (pp, new_seq) in enumerate(new_seqs):
        # s is tgt_seq, but expanded based on new_seq's dash pattern
        s: str = import_dash_pattern(already_dashed_seq=new_seq, 
            dest_seq_str=tgt_seq.sequence)
        seq_matching_to_pdb[new_seq.sequence] = simple_match(Sequence(s), new_seq, quiet=False)
  
    # Temp: we don't want only the best score, but maybe the top five.
    #best_score: int = max(seq_matching_to_pdb.values())
    #new_seqs = [(p, s) for (p, s) in new_seqs if s.sequence in seq_matching_to_pdb.keys() and seq_matching_to_pdb[s.sequence] == best_score]
    best_score: int = sorted(list(set(seq_matching_to_pdb.values())))[-20]
    new_seqs = [(p, s) for (p, s) in new_seqs if s.sequence in seq_matching_to_pdb.keys() and seq_matching_to_pdb[s.sequence] >= best_score]
    
    print("There are {n} sequences that match your template sequence to a score of {score}:".format(n=len(new_seqs), score=best_score))
    for ii, (p, q) in enumerate(new_seqs):
        print("\tSequence {serial} -- score {score}\n\t\tPDB: {pdb}\n\t\tSEQ: {seq}\n"
            .format(serial=ii, score=seq_matching_to_pdb[q.sequence], pdb=p, seq=q.sequence))

    return new_seqs

def match_pdb_to_best_sequence(pdb: str, MSA_file: str) -> List[Tuple[str, Sequence]]:
    """
    Get FASTA information from provided PDB. To which sequence is it closest? Pick it.
    """

    MSA_seqs = get_seqs_mfa(MSA_file)
    modomics_seq = modomics_from_pdb(pdb) # type Sequence
    print("Attempting to align", modomics_seq, "to MSA templates length", len(MSA_seqs))
    return match_seq_to_best_template_seq(modomics_from_pdb(pdb), MSA_seqs)

def add_dash_recursive(template: Sequence, trial: str, dashes, current_best):
    """
    Not in current use -- this is a very expensive function that 
    is good for aligning very difficult sequences. At the moment we
    have been using a few manual tweaks after automated alignment
    and that has been good enough.
    """
    
    def filled_trial_seq(trial: str, ii: int) -> str:
        return trial[:ii]+'-'+trial[ii:]+'-'*(len(template)-len(trial)-1)

    n_left_to_add = len(template)-len(trial)
    # Construct
    if len(dashes) == n_left_to_add:
        # We have enough.
        complete_trial_string = "-"*len(template)
        trial_index = 0
        for complete_index in range(len(template)-1):
            if complete_index in dashes: continue
            else:
                complete_trial_string = complete_trial_string[:complete_index]+trial[trial_index]+complete_trial_string[complete_index+1:]
                trial_index += 1
                if trial_index == len(trial): break
        score = simple_match(template, Sequence(complete_trial_string))
        if current_best is None or score > current_best[1]:
            current_best = (complete_trial_string, score)
    else:
        for ii in range(len(template)):
            if ii in dashes: continue
            else:
                new_dashes = list(dashes)
                new_dashes.append(ii)
            
                current_best = add_dash_recursive(template, trial, new_dashes, current_best)

    return current_best

def dash_positions(dashed_seq_str: str) -> List[int]:
    return [i for (i, c) in enumerate(dashed_seq_str) if c == '-']

def seqs_with_dashes(dashes: List[int], seq_str: str, max_len: int) -> List[str]:
    seqs = [seq_str]
    for dash in dashes:
        seq_str = seq_str[0:dash] + '-' + seq_str[dash:]
        if len(seq_str) > max_len: break
        seqs.append(seq_str)

    return seqs

def import_dash_pattern(already_dashed_seq: Sequence, dest_seq_str: str) -> str:
    """
    Put all the dashes from dashed_seq into other_seq.
    Try to be a LITTLE clever here. We don't want to just shove every in there
    in case there is an insertion. Maybe we should test each insertion to see
    if it improves alignment.

    Importantly, we have a maximum length to contend with...
    """

    dash_pos: List[int] = dash_positions(already_dashed_seq.sequence)
    possible_dest_seq_strs: List[str] = seqs_with_dashes(dash_pos, dest_seq_str, len(already_dashed_seq))
    score: int = None
    revised_seq_str: str = ""
    for possible_dest_seq_str in possible_dest_seq_strs:
        possible_dest_seq: Sequence = Sequence(possible_dest_seq_str)
        newscore: int = simple_match(possible_dest_seq, already_dashed_seq, quiet=True)
        if score is None or newscore > score:
            score = newscore
            revised_seq_str = possible_dest_seq_str  

    return revised_seq_str + (len(already_dashed_seq) - len(revised_seq_str)) * '-'

def align_template_library(MSA: str) -> None:
    """
    This was used to set up all the actual trna structure library sequences.
    """
    import glob, os
    f = open(my_loc() + "/data/all_trna_structure_seqs.dat", "w")

    for pdb in glob.glob(my_loc() + "/trnas/*.pdb"):
        print("Evaluating PDB:", pdb)
        seq = match_pdb_to_best_sequence(pdb, MSA)
        modomics_pdb_seq_str = modomics_from_pdb(pdb).sequence
        if seq is None:
            f.write("{}\t{} UNALIGNED_PDB_SEQ\n".format(pdb, modomics_pdb_seq_str))
            continue
        f.write("{}\t{} ALIGNED_TO\n".format(pdb, seq[0][1].sequence))
        f.write("{}\t{} PDB_SEQ\n".format(pdb, import_dash_pattern(seq[0][1], modomics_pdb_seq_str)))
        #exit()
        #print(pdb, seq)
        #print(pdb, import_dash_pattern(seq, modomics_from_pdb(pdb)))
    f.close()

def main(args):
    if args.pdb is not None: 
        pass
    
    import os
    if not os.path.exists(my_loc() + "/data/all_trna_structure_seqs.dat"):
        print("Regenerating aligned template library")
        align_template_library(my_loc() + "/data/all_trna.mfa")

    templates = get_seqs(my_loc() + "/data/all_trna_structure_seqs.dat")

    tgt_seq = Sequence("", 
        ".(((((((..((((...........)))).(((((.......)))))........................(((((.......))))))))))))....")
    with open(args.seq_file[0]) as f:
        tgt_seq.sequence = f.readlines()[0].strip()

    pdb_seq_list = match_seq_to_best_template_seq(tgt_seq, templates)
    
    # Maybe there are many returned! That's cool; do them all
    for p, s in pdb_seq_list:
        if 'a' in tgt_seq.sequence and 'g' in tgt_seq.sequence and 'c' in tgt_seq.sequence and 'u' in tgt_seq.sequence:
            # annotated seq format, must translate first.
            tgt_seq = ann_to_mod(tgt_seq)
        remodel_new_sequence(s, tgt_seq, p, args.nstruct, '', args.defer, args.aggressive)


"""
We need to make more stuff optional so that users can request that we use our pre-existing template library
"""

import argparse

if __name__ == '__main__':


    #print(repeat("simple_match(Sequence('ACUU--T-'), Sequence('CATP-GG-'))", number=10000, repeat=3, globals=globals()))
    
    parser = argparse.ArgumentParser(description="Thread tRNAs onto templates")
    parser.add_argument('--pdb', nargs='?', help='template (if omitted, use shipped library')
    parser.add_argument('--mapfile', nargs='?', help='electron density mapfile associated with template (useful for keeping minimization close)')
    parser.add_argument('--seq_file', nargs=1, help='target sequence file in modomics format')
    parser.add_argument('--nstruct', nargs='?', help='number of structures per template', default=1)
    parser.add_argument('--defer', nargs=1, help='defer execution (for cluster runs)', default=True)
    parser.add_argument('--aggressive', nargs=1, help='make aggressive deletions (every mismatch)', default=False)

    args = parser.parse_args()
    main(args)















