#!/usr/bin/env python3

import sys, subprocess

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
#   Data file manipulation: reading in sequence data and PDB file references.
###

def my_loc():
    import os
    return os.path.dirname(os.path.realpath(__file__))

def get_seqs_mfa(fn):
    """
    Get sequences out of the tRNAdb mfa
    """
    lines = []
    with open(fn) as f:
        lines = f.readlines()
    sequences = [l.strip() for l in lines[1::2]]
    print(len(sequences), "sequences")
    return sequences

def get_seqs(fn):
    lines = []
    with open(fn) as f:
        lines = f.readlines()
    sequences = {l.strip().split()[0]: l.strip().split()[1] for l in lines if "ALIGNED_TO" in l }
    #print(sequences)
    return sequences

###
#   Sequence translation (among Rosetta and Modomics nomenclatures)
###

# YG is yW, AET is MTA, G7M is 7MG, QUO is   Q, 1RN is CNS, YYG is yW [ish], 70U is MST, 12A is M26, 2MU is 52U, 6IA is I6A
tlc_to_mod = { #                                                                       note this is the enantiomer at the AA
               "GTP": 'G', " YG": 'Y', "AET": 'E', "G7M": '7', "QUO": 'Q', "1RN": '$', "YYG": 'Y', "  I": 'I',
               "70U": '3', "12A": '[', "2MU": '\\', "6IA": '+',

               "  A": 'A', "  C": 'C', "  G": 'G', "  U": 'U', "1MA": '"', "2MA": '/', "I6A": '+', "MIA": '*',
               "MA6": '=', "T6A": '6', "MTA": 'E', "M26": '[', "OMA": ':', "INO": 'I', "M1I": 'O', "RIA": '^',
               "IHA": '`', "S2C": '%', "OMC": 'B', "A4C": 'M', "5MC": '?', "3MC": '\'', "K2C": '}', "F5C": '>',
               "52C": '<', "1MG": 'K', "2MG": 'L', "OMG": '#', "M2G": 'R', "M3G": '|', "7MG": '7', "FA7": '(',
               "  Q": 'Q', " yW": 'Y', "O2W": 'W', "2SU": '2', "OMU": 'J', "4SU": '4', "5MU": 'T', "2ST": 'F',
               "52U": '\\', "NMT": '{', "NST": 'S', "CST": '&', "CMT": '~', "MOT": '1', "MST": '3', "OAU": 'V',
               "MOU": '5', "CNT": '!', "CNS": '$', "CNM": ')', "APU": 'X', "MHU": ',', "H2U": 'D', "PSU": 'P',
               "1PU": ']', "MPU": 'Z' }
mod_to_tlc = dict(zip(tlc_to_mod.values(), tlc_to_mod.keys()))

def mod_from_tlc(tlc):
    """
    Translates modomics single letter codes from 3LCs
    """
    if tlc in tlc_to_mod:
        return tlc_to_mod[tlc]

    print("Unrecognized", tlc)
    exit()

def modomics_from_pdb(pdb):
    pdblines = []
    modomics_seq = ""
    with open(pdb) as f:
        pdblines = f.readlines()
    for l in pdblines:
        if " C4'" in l: modomics_seq += mod_from_tlc(l[17:20])
    return modomics_seq

def ann_to_mod(ann_seq):
    mod_seq = ""
    for c in ann_seq:
        if c == 'a': mod_seq += 'A'
        elif c == 'c': mod_seq += 'C'
        elif c == 'g': mod_seq += 'G'
        elif c == 'u': mod_seq += 'U'
        else:
            assert(c == 'X')
            print("not yet supported")
            # AMW TODO
            exit()

    return mod_seq

###
#   Sequence alignment scoring
###

def simple_match(seq1, seq2):
    """
    Assumes seq1 and seq2 are aligned already (i.e., a 'static score').
    Rewards similarity in length, letter matches, purine matches, and U variant
    matches. Could reward other variants! Those just took priority because of 
    the ultra-characteristic TP sequence.

    Ooh, consider not penalizing complementary WC mismatches within helices?
    would need SS
    """

    score = len(seq1) - len(seq2)
    for c1, c2 in zip(seq1, seq2):
        if c1 == c2: score += 5
        else:
            if c1 == 'A' and c2 == 'G': score += 2
            elif c1 == 'G' and c2 == 'A': score += 2
            elif c1 == 'C' and c2 == 'U': score += 2
            elif c1 == 'U' and c2 == 'C': score += 2
            elif c1 == 'U' and c2 == 'T': score += 3
            elif c1 == 'T' and c2 == 'U': score += 3
            elif c1 == 'U' and c2 == 'P': score += 3
            elif c1 == 'P' and c2 == 'U': score += 3
            elif (c1 == '-' or c2 == '-'): score -= 10

    #print(score, seq1, seq2)
    return score


"""
DOMAIN ANALYSIS:

1. Proximal stuff
2. Never rebuild key tert contact (GG t-loop) unless it is directly altered 
"""

# This string 'colors' bits of tRNA structure that should be rebuilt together.
trna_structure_coloring = "P(((((((VV((((DDDDD..DDDD))))V(((((SSSSSSS)))))VVVVVVVVVVVVVVVVVVVVVVVV(((((.....DD))))))))))))PPPP"

def match_seq_to_best_template_seq(tgt_seq, templates):
    
    print("There are in total", len(templates))
    # Filter sequences for those of same true length. If 1, done.
    # If zero, quit entirely; no point.
    new_seqs = [(p,s) for (p,s) in templates.items()]
    ii = 0
    #while len(new_seqs) == 0:
    #    new_seqs = [(p,s) for (p,s) in templates.items() if abs(len(tgt_seq.replace('-', '')) - len(s.replace('-', ''))) < ii]
    #    ii += 1
    #for _, seq in new_seqs:
    #    print(seq)
    #    print(tgt_seq)
    #print("There are length-matched", len(new_seqs))
    #if len(new_seqs) == 0: return None#exit(0)

    pdb_len = len(new_seqs[0][1].replace('-',''))

    # Sort sequences by number of matching characters (exact + purine/pyr+TP/U),
    # take best count only. If 1, done.
    seq_matching_to_pdb = { s: simple_match(tgt_seq, s.replace('-','')) for _, s in new_seqs }
    best_score = max(seq_matching_to_pdb.values())
    new_seqs = [(p, s) for (p, s) in new_seqs if seq_matching_to_pdb[s] == best_score]
    print("There are, for top seq exact-match", len(new_seqs), "(each scores {score} out of {length})".format(score=best_score, length=pdb_len))
    for p,s in new_seqs:
        print(p)
        print(s)

    return new_seqs

def match_pdb_to_best_sequence(pdb, MSA_file):
    """
    Get FASTA information from provided PDB. To which sequence is it closest? Pick it.
    """

    MSA_seqs = get_seqs_mfa(MSA_file)
    # Has to have some dummy key value; must be unique.
    templates = {s: s for s in MSA_seqs}
    modomics_seq = modomics_from_pdb(pdb)
    print("Attempting to align", modomics_seq, "to MSA templates length", len(templates))
    return match_seq_to_best_template_seq(modomics_from_pdb(pdb), templates)




def thread_sequence_on_pdb(seq, pdb, mapfile):
    """
    Threads a sequence using rna_thread_and_minimize. Could use pyrosetta for this, I imagine.
    Currently uses an un-translated modomics sequence; could imagine using annotated sequence
    (if we translate ourselves).


    This was appropriate back when we were INITIALLY doing some threading just to get into the
    language of the old MFA. That is dumb and unnecessary.
    """

    # If we already have a sequence match, don't bother threading. What really matters is our careful sequence
    # matching allowed us to align to the overall alignment from the alignment file.
    print("PYTHON: threading {seq} onto {pdb_seq}.".format(seq=seq.replace('-', ''), pdb_seq=modomics_from_pdb(pdb)))
    if seq.replace('-', '') == modomics_from_pdb(pdb):
        # Copy input PDB to where it would be expected as output.
        subprocess.run(['cp', pdb, pdb.replace('.pdb', '_0001.pdb')])
        return pdb.replace('.pdb', '_0001.pdb')

    # Deal with one issue: reading an absolute path template library, copy to PWD.
    if "/" in pdb:
        # Copy input PDB to where it would be expected as output.
        subprocess.run(['cp', pdb, './template.pdb'])
        pdb = 'template.pdb'

    # Question of how to incorporate native constraints, a guaranteed foldtree, cart bonded, density, etc.
    # unresolved.
    command = ["/home/andy/Rosetta/main/source/bin/rna_thread_and_minimize", 
        "-s", pdb, 
        "-seq", seq.replace('-',''), 
        "-input_sequence_type", "MODOMICS", 
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
        "-ignore_zero_occupancy", "false",
        "-guarantee_no_DNA", "true",
        "-set_weights", "other_pose", "0.0", "intermol", "0.0", "loop_close", "0.0", "free_suite", "0.0", "free_2HOprime", "0.0"]
    if mapfile is not None:
        command.append("elec_dens_fast")
        command.append("10.0")
        command.append("-edensity:mapfile")
        command.append(mapfile)
    print(command)
    subprocess.run(command)

    # Return new name -- helps ensure that we keep track properly.
    return pdb.replace('.pdb', '_0001.pdb')






def remodel_new_sequence(seq, tgt_seq, pdb, mapfile=None):
    """
    This is the main workhorse. We need to figure out a path from A to B. How do we do this?
    1. Calculate common structure. This is more restrictive than 'common residues' because we have to eliminate whole loops
    or loop segments after mutation.
    2. Trim input PDB to common structure. Figure out which sequence positions this is.
    3. Thread new sequence using above effort onto common structure. [How to restrict minimization?]
    4. 
    """
    
    def remove_if_of_common_color(mm, common_structure, trna_structure_coloring=trna_structure_coloring):
        common_color_set = [d for d in range(len(trna_structure_coloring)) if trna_structure_coloring[d] == trna_structure_coloring[mm]]
        for c in common_color_set: common_structure.remove(c)
        return common_structure

    print(seq)
    print(tgt_seq)
    
    # Deal with one issue: reading an absolute path template library, copy to PWD.
    old_pdb = pdb
    if "/" in pdb:
        # Copy input PDB to where it would be expected as output.
        subprocess.run(['cp', pdb, './template.pdb'])
        old_pdb = pdb[-10:-4]
        pdb = 'template.pdb'

    mismatches = []
    for i, (c1, c2) in enumerate(zip(seq, tgt_seq)):
        if c1 != c2 and ((c1 == '-' and c2 != '-') or (c2 == '-' and c1 != '-')):
            print("mismatch at {pos}".format(pos=i))
            mismatches.append(i)
    common_structure = [i for i in range(len(seq))]
    for mm in mismatches:
        if mm in common_structure:
            common_structure = remove_if_of_common_color(mm, common_structure)
    print("So, common structure remaining is", common_structure)


    # Trim input PDB containing only common_structure residues
    pdb_pos2align_pos = {}
    pdb_pos = 1
    for align_pos in range(len(seq)):
        if seq[align_pos] == '-': continue
        pdb_pos2align_pos[pdb_pos] = align_pos
        pdb_pos += 1
    align_pos2pdb_pos = dict(zip(pdb_pos2align_pos.values(), pdb_pos2align_pos.keys()))

    def trim_pdb_to_res(pdb, new_pdb, pos):
        def nres(pdb):
            with open(pdb) as f:
                return len([ 0 for l in f.readlines() if ' P  ' in l])

        def atom(line):
            return line[12:16]
        def res(line):
            return line[22:26]
        lines = []
        with open(pdb) as f:
            lines = f.readlines()
        # Renumber the PDB, get a dict from seqpos => lines
        # write all in pos
        subprocess.run(['cp', pdb, new_pdb])
        print("I think there are {} res".format(nres(new_pdb)))
        subprocess.run(['renumber_pdb_in_place.py', new_pdb, "A:1-{}".format(nres(new_pdb))])
        with open(new_pdb) as f:
            lines = f.readlines()
        res_to_line_dict = { x: [l for l in lines if int(res(l).strip()) == x and int(res(l).strip()) in pos] for x in pos }
        with open(new_pdb, "w") as f:
            for x in pos:
                for l in res_to_line_dict[x]:
                    f.write(l)
            #res_index = 1
            #lineres = None
            #for line in lines:
            #    if line[0:4] != "ATOM" and line[0:6] != "HETATM": continue
            #    #print(atom(line))
            #    #print(line, res_index)
            #    if res_index in pos: f.write(line)
            #    if res(line) != lineres:
            #        if lineres is not None: res_index += 1
            #        lineres = res(line)

    pdb_pos_trim = [align_pos2pdb_pos[c] for c in common_structure if c in align_pos2pdb_pos.keys()]
    print(pdb_pos_trim)
    trim_pdb_to_res(pdb, pdb.replace('.pdb', '_trimmed.pdb'), pdb_pos_trim)

    # Thread them to their new identities.
    tgt_seq_for_thread = "".join([tgt_seq[c] for c in common_structure])

    command = ["/home/andy/Rosetta/main/source/bin/rna_thread_and_minimize", 
        "-s", pdb.replace('.pdb', '_trimmed.pdb'), 
        "-seq", tgt_seq_for_thread.replace('-',''), 
        "-input_sequence_type", "MODOMICS", 
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
        "-guarantee_no_DNA", "true",
        "-ignore_zero_occupancy", "false",
        "-overwrite",
        "-set_weights", "other_pose", "0.0", "intermol", "0.0", "loop_close", "0.0", "free_suite", "0.0", "free_2HOprime", "0.0"]
    # Testing not using mapfile for the new no-bb motion option
    #if mapfile is not None:
    #    command.append("elec_dens_fast")
    #    command.append("10.0")
    #    command.append("-edensity:mapfile")
    #    command.append(mapfile)
    print(command)
    subprocess.run(command)
    
    #subprocess.run(["mv", pdb.replace('.pdb', '_trimmed_0001.pdb'), "input.pdb"])
    # Remove LINK records. They will be misleading because of a bug in the threading code.
    with open(pdb.replace('.pdb', '_trimmed_0001.pdb')) as infile:
        with open("{}_input.pdb".format(old_pdb), "w") as outfile:
            for line in infile.readlines():
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    outfile.write(line)

    old_seq2new_seq = {}
    oi, ni = 0, 0
    for oc, nc in zip(seq, tgt_seq):
        if oc != '-': oi += 1
        if nc != '-': ni += 1
        
        old_seq2new_seq[oi] = ni

    # AMW: incidentally we are starting from T:1-76
    # We would have to check that this is robust to (say) T:10-85
    new_trimmed = [old_seq2new_seq[p] for p in pdb_pos_trim]

    def dashed(positions, chain):
        #print(positions)
        acc_start = None
        acc_end = None
        numbering = ""
        for p in range(1, max(positions)+1):
            #print("examining", p)
            if p in positions and p != max(positions):
                #print("\tfound in positions", p)
                if acc_start is None:
                    acc_start = p
                else:
                    acc_end = p
                #print("now acc_start acc_end", acc_start, acc_end)
            elif p != max(positions):
                if acc_start is not None and acc_end is not None:
                    numbering += "{ch}:{start}-{end} ".format(ch=chain, start=acc_start, end=acc_end)
                    acc_start = None
                    acc_end = None
            else:
                if p in positions:
                    acc_end = p
                numbering += "{ch}:{start}-{end} ".format(ch=chain, start=acc_start, end=acc_end)

        return numbering

    numbering = dashed(new_trimmed, 'A')
    print(numbering)
    
    # renumber input if needed? For sure to new chain.
    # ugly -- but numbering has spaces in it
    command = ["renumber_pdb_in_place.py", "{}_input.pdb".format(old_pdb)]
    command.extend(numbering.split())
    print(command)
    subprocess.run(command)

    # Create the overall fasta.

    def annotated_seq_of(mod_seq):
        ann_seq = ""
        for c in mod_seq:
            if c == 'A': ann_seq += 'a'
            elif c == 'C': ann_seq += 'c'
            elif c == 'G': ann_seq += 'g'
            elif c == 'U': ann_seq += 'u'
            else: ann_seq += 'X[' + mod_to_tlc[c] + ']'
        return ann_seq

    with open("target.fasta", "w") as f:
        print(tgt_seq)
        print(tgt_seq.replace('-',''))
        print(len(tgt_seq.replace('-','')))
        f.write(">foo A:1-{}\n".format(len(tgt_seq.replace('-', ''))))
        f.write("{}\n".format(annotated_seq_of(tgt_seq.replace('-', ''))))


    # Do denovo run.
    command = ["/Users/amw579/dev_ui/Rosetta/main/source/bin/rna_denovo", 
        "-s", "{}_input.pdb".format(old_pdb), 
        "-fasta", "target.fasta", 
        "-minimize_rna", "true",
        "-include_neighbor_base_stacks", "true",
        "-motif_mode", "true",
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
        "-set_weights", "other_pose", "0.0", "intermol", "0.0", "loop_close", "0.0", "free_suite", "0.0", "free_2HOprime", "0.0",
        "-guarantee_no_DNA", "true",
        "-use_legacy_job_distributor", "true",
        "-nstruct", "1",
        "-out:file:silent", "{}_based_modeled.out".format(old_pdb)]
    
    subprocess.run(command)
    subprocess.run(["extract_lowscore_decoys.py", "{}_based_modeled.out".format(old_pdb)])



def import_dash_pattern(dashed_seq, other_seq):
    """
    Put all the dashes from dashed_seq into other_seq.
    Try to be a LITTLE clever here. We don't want to just shove every in there
    in case there is an insertion. Maybe we should test each insertion to see
    if it improves alignment.
    """
    def dash_positions(dashed_seq):
        return [i for (i, c) in enumerate(dashed_seq) if c == '-']

    def seqs_with_dashes(dashes, seq):
        seqs = [seq]
        for i in range(len(dashes)):
            seq = seq[0:dashes[i]] + '-' + seq[dashes[i]:]
            seqs.append(seq)

        return seqs

    dash_pos = dash_positions(dashed_seq)
    seqs = seqs_with_dashes(dash_pos, other_seq)
    score = None
    for seq in seqs:
        newscore = simple_match(seq, dashed_seq)
        #print(newscore, seq)
        if score is None or newscore > score:
            score = newscore
            other_seq = seq  

    #for dash in dash_pos:
    #    proposed_new_seq = other_seq[0:dash] + '-' + other_seq[dash:]
    #    match_scores = match(proposed_new_seq, dashed_seq), match(other_seq, dashed_seq)
    #    if match_scores[0] >= match_scores[1]:
    #        other_seq = proposed_new_seq

    #print(other_seq)

    return other_seq + (len(dashed_seq) - len(other_seq)) * '-'

def align_template_library(MSA):
    """
    This was used to set up all the actual trna structure library sequences.
    """
    import glob, os
    f = open(my_loc() + "/data/all_trna_structure_seqs.dat", "w")

    for pdb in glob.glob(my_loc() + "/trnas/*.pdb"):
        print(pdb)
        seq = match_pdb_to_best_sequence(pdb, MSA)
        if seq is None:
            f.write("{}\t{} UNALIGNED_PDB_SEQ\n".format(pdb, modomics_from_pdb(pdb)))
            continue
        f.write("{}\t{} ALIGNED_TO\n".format(pdb, seq[0][1]))
        f.write("{}\t{} PDB_SEQ\n".format(pdb, import_dash_pattern(seq[0][1], modomics_from_pdb(pdb))))

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

    with open(args.seq_file[0]) as f:
        tgt_seq = f.readlines()[0].strip()

    pdb_seq_list = match_seq_to_best_template_seq(tgt_seq, templates)
    
    # Maybe there are many returned! That's cool; do them all
    for p, s in pdb_seq_list:
        if 'a' in tgt_seq and 'g' in tgt_seq and 'c' in tgt_seq and 'u' in tgt_seq:
            # annotated seq format, must translate first.
            tgt_seq = ann_to_mod(tgt_seq)
        remodel_new_sequence(s, tgt_seq, p, args.nstruct)


"""
We need to make more stuff optional so that users can request that we use our pre-existing template library
"""

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Thread tRNAs onto templates")
    parser.add_argument('--pdb', nargs='?', help='template (if omitted, use shipped library')
    parser.add_argument('--mapfile', nargs='?', help='electron density mapfile associated with template (useful for keeping minimization close)')
    parser.add_argument('--seq_file', nargs=1, help='target sequence file in modomics format')
    parser.add_argument('--nstruct', nargs='?', help='number of structures per template', default=1)

    args = parser.parse_args()
    main(args)















