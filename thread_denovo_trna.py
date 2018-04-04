#!/usr/bin/env python3

import sys, subprocess

"""
OK, suppose that you have a file of the form:

>tdbR00000589|Aeropyrum_pernix|56636|Arg|UCU
-GGGCCCGUAGCUCAGCCAGGAC-AGAGCGCCGGCCUUCUAAGCCGGUG-------------------CUGCCGGGUUCAAAUCCCGGCGGGCCCGCCA
>tdbR00000588|Aeropyrum_pernix|56636|Asp|GUC
-GCCGCGGUAGUAUAGCCUGGACUAGUAUGCGGGCCUGUCAAGCCCGUG-------------------A-CCCGGGUUCAAAUCCCGGCCGCGGCGCCA
>tdbR00000583|Aeropyrum_pernix|56636|Cys|GCA
-GCCGGGGUGGCCGAGC--GGUCUAAGGCGGCGGGCUGCAGACCCGUUA-------G-----------UUCCCGGGUUCGAAUCCCGGCCCCGGCUCCA
 (((((((..((((...........)))).(((((.......)))))........................(((((.......))))))))))))....

i.e., a tRNA MSA. A plausible task is:

"I have a tRNA scaffold that looks like entry #2 and I need to turn it into a PDB that looks like 
entry #3." So we break this into:
    1. threading a (possibly NCNT containing) sequence onto a MATCHING scaffold (use density if
    available, cartesian minimization for sure). (to ensure that PDB truly matches entry #2)
    2. Use above strings w/ 'known secstruct' to figure out a 'domain decomposition' for a denovo
    job.
    3. Use above strings w/ knowledge about base pairing (i.e. -- eliminate everything but acgu)
    to create new secstruct string.


DOMAIN ANALYSIS:

If a modification (ins/del) happens in a position marked below with S you just need a simple loop
rebuild, b/c no tertiary contact

(((((((..((((...........)))).(((((SSSSSSS)))))........................(((((.......))))))))))))....

But those are unlikely. More likely are D-loop:

(((((((..((((DDDDDDDDDDD)))).(((((.......)))))........................(((((.......))))))))))))....

in which case rebuild every residue EXCEPT G18-19 but add in 59-60 (between bp 53-61 and UA handle)

(((((((..((((DDDDD..DDDD)))).(((((.......)))))........................(((((.....DD))))))))))))....


or V-region (44-48, but can get insertions making everything longer ofc)

(((((((..((((...........)))).(((((.......)))))VVVVVVVVVVVVVVVVVVVVVVVV(((((.......))))))))))))....

makes some contacts:

(((((((VV((((...........))))V(((((.......)))))VVVVVVVVVVVVVVVVVVVVVVVV(((((.......))))))))))))....

Let's start with that assumption, and plot out what has to happen based on "- becomes letter" and
reverse.
"""


# This string 'colors' bits of tRNA structure that should be rebuilt together.
trna_structure_coloring = "-(((((((VV((((DDDDD..DDDD))))V(((((.......)))))VVVVVVVVVVVVVVVVVVVVVVVV(((((.....DD))))))))))))...."


def get_seqs(fn):
    lines = []
    with open(fn) as f:
        lines = f.readlines()
    sequences = [l.strip() for l in lines[1::2]]
    #print(sequences)
    return sequences

tlc_to_mod = { "  A": 'A', "  C": 'C', "  G": 'G', "  U": 'U', "1MA": '"', "2MA": '/', "I6A": '+', "MIA": '*',
               "MA6": '=', "T6A": '6', "MTA": 'E', "M26": '[', "OMA": ':', "INO": 'I', "M1I": 'O', "R1A": '^',
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

    print("Unrecognized")
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

    return mod_seq


def match_pdb_to_best_sequence(pdb, sequences):
    """
    Get FASTA information from provided PDB. To which sequence is it closest? Pick it.
    """

    def nres_pdb(pdb):
        pdblines = []
        i = 0
        with open(pdb) as f:
            pdblines = f.readlines()
        for l in pdblines:
            if " C4'" in l: i += 1
        return i

    def len_match(pdb, s):
        return len(s.replace('-', '')) == nres_pdb(pdb)

    def n_match_char(s1, s2):
        i = 0
        for c1, c2 in zip(s1, s2):
            if c1 == c2: i += 1
        return i 

    print("There are in total", len(sequences))
    # Filter sequences for those of same true length. If 1, done.
    # If zero, quit entirely; no point.
    sequences = [s for s in sequences if len_match(pdb, s.replace('-',''))]
    print("There are length-matched", len(sequences))
    if len(sequences) == 0: exit(0)

    pdb_len = len(sequences[0].replace('-',''))

    pdb_seq = modomics_from_pdb(pdb)
    # Sort sequences by number of matching characters (exact),
    # take best count only. If 1, done.
    seq_matching_to_pdb = { s: n_match_char(pdb_seq, s.replace('-','')) for s in sequences }
    best_score = max(seq_matching_to_pdb.values())
    sequences = [s for s in sequences if seq_matching_to_pdb[s] == best_score]
    print("There are, for top seq exact-match", len(sequences), "(each scores {score} out of {length})".format(score=best_score, length=pdb_len))

    if best_score == pdb_len:
        # All are perfect; no need to distinguish.
        return sequences[0]

    # Purine/pyrimidine sequence matching.
    #for seq in sequences:
    #    print(len(seq.replace('-','')), seq)


    return sequences[0]

def thread_sequence_on_pdb(seq, pdb, mapfile):
    """
    Threads a sequence using rna_thread_and_minimize. Could use pyrosetta for this, I imagine.
    Currently uses an un-translated modomics sequence; could imagine using annotated sequence
    (if we translate ourselves).
    """

    # If we already have a sequence match, don't bother threading. What really matters is our careful sequence
    # matching allowed us to align to the overall alignment from the alignment file.
    print("PYTHON: threading {seq} onto {pdb_seq}.".format(seq=seq.replace('-', ''), pdb_seq=modomics_from_pdb(pdb)))
    if seq.replace('-', '') == modomics_from_pdb(pdb):
        # Copy input PDB to where it would be expected as output.
        subprocess.run(['cp', pdb, pdb.replace('.pdb', '_0001.pdb')])
        return


    # Question of how to incorporate native constraints, a guaranteed foldtree, cart bonded, density, etc.
    # unresolved.
    command = ["/Users/amw579/dev_ui/Rosetta/main/source/bin/rna_thread_and_minimize", 
        "-s", pdb, 
        "-seq", seq.replace('-',''), 
        "-input_sequence_type", "MODOMICS", 
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
        "-set_weights", "other_pose", "0.0", "intermol", "0.0", "loop_close", "0.0", "free_suite", "0.0", "free_2HOprime", "0.0"]
    if mapfile is not None:
        command.append("elec_dens_fast")
        command.append("10.0")
        command.append("-edensity:mapfile")
        command.append(mapfile)
    print(command)
    subprocess.run(command)

def remodel_new_sequence(seq, tgt_seq, pdb, mapfile):
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
        def atom(line):
            return line[12:16]
        lines = []
        with open(pdb) as f:
            lines = f.readlines()
        with open(new_pdb, "w") as f:
            res_index = 1
            for line in lines:
                if line[0:4] != "ATOM" and line[0:6] != "HETATM": continue
                #print(atom(line))
                #print(line, res_index)
                if res_index in pos: f.write(line)
                if atom(line) == " P  ": res_index += 1

    pdb_pos_trim = [align_pos2pdb_pos[c] for c in common_structure if c in align_pos2pdb_pos.keys()]
    print(pdb_pos_trim)
    trim_pdb_to_res(pdb, pdb.replace('.pdb', '_trimmed.pdb'), pdb_pos_trim)

    # Thread them to their new identities.
    tgt_seq_for_thread = "".join([tgt_seq[c] for c in common_structure])

    command = ["/Users/amw579/dev_ui/Rosetta/main/source/bin/rna_thread_and_minimize", 
        "-s", pdb.replace('.pdb', '_trimmed.pdb'), 
        "-seq", tgt_seq_for_thread.replace('-',''), 
        "-input_sequence_type", "MODOMICS", 
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
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
        with open("input.pdb", "w") as outfile:
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
    command = ["renumber_pdb_in_place.py", "input.pdb"]
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
        f.write(">foo A:1-{}\n".format(len(tgt_seq.replace('-', ''))))
        f.write(annotated_seq_of(tgt_seq.replace('-', '')))


    # Do denovo run.
    command = ["/Users/amw579/dev_ui/Rosetta/main/source/bin/rna_denovo", 
        "-s", "input.pdb", 
        "-fasta", "target.fasta", 
        "-minimize_rna", "true",
        "-include_neighbor_base_stacks", "true",
        "-motif_mode", "true",
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
        "-set_weights", "other_pose", "0.0", "intermol", "0.0", "loop_close", "0.0", "free_suite", "0.0", "free_2HOprime", "0.0",
        "-nstruct", "100",
        "-out:file:silent", "modeled.out"]
    
    subprocess.run(command)
    subprocess.run(["extract_lowscore_decoys.py", "modeled.out"])





def main(pdb="start.pdb", mapfile=None, tgt_seq=None):
    import os
    dirname = os.path.dirname(__file__)
    sequences = get_seqs(dirmame+"data/all_trna.mfa")
    seq = match_pdb_to_best_sequence(pdb, sequences)
    # There are no gaps with this sequence.
    thread_sequence_on_pdb(seq, pdb, mapfile)

    # Whatever the new sequence is, it's auto-aligned to seq and so we know how to set up a denovo job for it.
    if tgt_seq is not None:
        if 'a' in tgt_seq and 'g' in tgt_seq and 'c' in tgt_seq and 'u' in tgt_seq:
            # annotated seq format, must translate first.
            tgt_seq = ann_to_mod(tgt_seq)
        remodel_new_sequence(seq, tgt_seq, pdb.replace('.pdb', '_0001.pdb'), mapfile)

if __name__ == '__main__':
    main(*sys.argv[1:])
