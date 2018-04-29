import mypy
from typing import List, Tuple, Dict
import subprocess
from lib.Sequence import Sequence
from lib.str_manip import modomics_from_pdb, all_bp_partners, dashed, annotated_seq_of
from lib.file_io import exe, my_loc

def just_move(seq, pdb):
    if seq.replace('-', '') == modomics_from_pdb(pdb):
        # Copy input PDB to where it would be expected as output.
        subprocess.run(['cp', pdb, pdb.replace('.pdb', '_0001.pdb')])
        return True
    else: return False

def correct_for_template_path(pdb):
    # Deal with one issue: reading an absolute path template library, copy to PWD.
    if "/" in pdb:
        # Copy input PDB to where it would be expected as output.
        subprocess.run(['cp', pdb, './template.pdb'])
        pdb = 'template.pdb'
    return pdb





def atom(line: str) -> str:
    return line[12:16]
def res(line: str) -> str:
    return line[22:26]
 
def trim_pdb_to_res(pdb: str, new_pdb: str, pos: List[int]) -> None:
    def nres(pdb: str) -> int:
        n = 0
        lastnum = None
        with open(pdb) as f:
            for l in f.readlines():
                if res(l) != lastnum:
                    lastnum = res(l)
                    n += 1
            #return len([ 0 for l in f.readlines() if ' P  ' in l])
        return n

    lines = [] # type: List[str]
    with open(pdb) as f:
        lines = f.readlines()
    # Renumber the PDB, get a dict from seqpos => lines
    # write all in pos
    subprocess.run(['cp', pdb, new_pdb])
    print("I think there are {} res".format(nres(new_pdb)))
    subprocess.run(['renumber_pdb_in_place.py', new_pdb, "A:1-{}".format(nres(new_pdb))])
    #exit()
    with open(new_pdb) as f:
        lines = f.readlines()
    res_to_line_dict = { x: [l for l in lines if int(res(l).strip()) == x and int(res(l).strip()) in pos] for x in pos }
    #for k in res_to_line_dict.keys():
    #    for l in res_to_line_dict[k]:
    #        print("{}: {}".format(k, l.strip()))
    #exit()
    with open(new_pdb, "w") as f:
        for x in pos:
            for l in res_to_line_dict[x]:
                f.write(l)

def cull_LINKs(in_pdb: str, out_pdb: str) -> None:
    with open(in_pdb) as infile:
        with open(out_pdb, "w") as outfile:
            for line in infile.readlines():
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    outfile.write(line)





def thread_sequence_on_pdb(seq: str, pdb: str, mapfile: str) -> str:
    """
    Threads a sequence using rna_thread_and_minimize. Could use pyrosetta for this, I imagine.
    Currently uses an un-translated modomics sequence; could imagine using annotated sequence
    (if we translate ourselves).


    This was appropriate back when we were INITIALLY doing some threading just to get into the
    language of the old MFA. That is dumb and unnecessary.
    """

    # If we already have a sequence match, don't bother threading. What really matters is our careful sequence
    # matching allowed us to align to the overall alignment from the alignment file.
    print("PYTHON: threading {seq} onto {pdb_seq}.".format(seq=seq.replace('-', ''), pdb_seq=modomics_from_pdb(pdb).sequence))
    
    if just_move(seq, pdb):
        return pdb.replace('.pdb', '_0001.pdb')

    pdb = correct_for_template_path(pdb)

    # Question of how to incorporate native constraints, a guaranteed foldtree, cart bonded, density, etc.
    # unresolved.
    command = [exe("rna_thread_and_minimize"),
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


def seq_map(seq1: str, seq2: str) -> Dict[int, int]:
    o2n: Dict[int, int] = {}
    oi, ni = 0, 0
    for oc, nc in zip(seq1, seq2):
        if oc != '-': oi += 1
        if nc != '-': ni += 1
        
        o2n[oi] = ni
    return o2n

def translation_dicts(seq_str: str, pdb_seq_str: str) -> Tuple[Dict[int, int], Dict[int, int]]:
    pdb_pos2align_pos = {}
    pdb_pos = 1
    try:
        for align_pos, align_char in enumerate(seq_str):
            if align_char == '-': continue
            if align_char != pdb_seq_str[pdb_pos-1]: 
                #pdb_pos += 1
                continue
            pdb_pos2align_pos[pdb_pos] = align_pos
            pdb_pos += 1
    except:
        print("ERROR:")
        print("PDB seq:", pdb_seq_str)
        print("SEQ seq:", seq_str)
    align_pos2pdb_pos = dict(zip(pdb_pos2align_pos.values(), pdb_pos2align_pos.keys()))
    return (pdb_pos2align_pos, align_pos2pdb_pos)

# This string 'colors' bits of tRNA structure that should be rebuilt together.
trna_structure_coloring = "P(((((((VV((((DDDDD..DDDD))))V(((((SSSSSSS)))))VVVVVVVVVVVVVVVVVVVVVVVV(((((.....DD))))))))))))PPPP"

# TODO: tests.
def remove_if_of_common_color(mm: int, common_structure: List[int], trna_structure_coloring=trna_structure_coloring) -> List[int]:
    """
    If one residue is deleted or added, what should be 
    rebuilt? Use the 'common coloring' string above to
    figure out what. Common letters -- as well as base 
    pair partners of modified residues -- get removed.
    """

    all_bp: Dict[int, int] = all_bp_partners(trna_structure_coloring)

    bp_chars = ['(', ')', '[', ']', '{', '}']
    positions_of_same_color: List[int] = []
    char_to_match = trna_structure_coloring[mm]
    if char_to_match == '.': 
        positions_of_same_color = [mm]
    elif char_to_match in bp_chars:
        positions_of_same_color = [mm, all_bp[mm]]
    else: positions_of_same_color = [ii for ii, cc in enumerate(trna_structure_coloring) 
            if cc == char_to_match]
    
    for cc in positions_of_same_color: 
        if cc in common_structure: common_structure.remove(cc)
    
    return common_structure


def get_common_structure(seq: Sequence, tgt_seq: Sequence, aggro: bool) -> List[int]:
    """
    Identify additions and deletions relative to the alignment,
    then separately account for the 'common structure' requirement
    so that we rebuild groups of interacting residues.
    """

    mismatches = []
    for ii, (c1, c2) in enumerate(zip(seq.sequence, tgt_seq.sequence)):
        if c1 != c2 and ((aggro and seq.secstruct[ii] == '.') or ((c1 == '-' and c2 != '-') or (c2 == '-' and c1 != '-'))):
            print("mismatch at {pos}".format(pos=ii))
            mismatches.append( ii )
    common_structure = [ii for ii, _ in enumerate(seq.sequence)]
    for mm in mismatches:
        if mm in common_structure:
            common_structure = remove_if_of_common_color(mm, common_structure)
    return common_structure


    
def write_fasta(fasta_file: str, tgt_seq: Sequence) -> None:    
    with open(fasta_file, "w") as f:
        #print(tgt_seq)
        #print(tgt_seq.sequence.replace('-',''))
        #print(len(tgt_seq.sequence.replace('-','')))
        f.write(">foo A:1-{}\n".format(len(tgt_seq.sequence.replace('-', ''))))
        f.write("{}\n".format(annotated_seq_of(tgt_seq.sequence.replace('-', ''))))


def remodel_new_sequence(seq: Sequence, tgt_seq: Sequence, pdb: str, nstruct: int, mapfile: str, defer: bool, aggro: bool) -> None:
    """
    This is the main workhorse. We need to figure out a path from A to B. How do
    we do this?
    1. Calculate common structure. This is more restrictive than 'common 
    residues' because we have to eliminate whole loops or loop segments after
    mutation.
    2. Trim input PDB to common structure. Figure out which sequence positions 
    this is.
    3. Thread new sequence using above effort onto common structure. [How to 
    restrict minimization? Let's say we just do chis for now, unless a mapfile
    is provided.]
    4. FARFAR or SWM run to build new residues.
    """

    print("PDB:", pdb)
    print("\nModeling from seq to target sequence:")
    print("SEQ:\n", seq)
    print("TGT:\n", tgt_seq)
    
    # Deal with one issue: reading an absolute path template library, copy to PWD.
    old_pdb: str = pdb[-10:-4]
    pdb = correct_for_template_path(pdb)

    # Now relies on the whole Sequence because we need to know if a position is a
    # loop (where, if aggro, a mere letter mismatch requires rebuild)
    common_structure = get_common_structure(seq, tgt_seq, aggro)
    
    print("So, common structure remaining (in alignment numbering) is", common_structure)
    print("COL:", trna_structure_coloring)
    # Sequence given common structure remaining
    print("SEQ:", "".join([c if i in common_structure else '-' for i,c in enumerate(seq.sequence) ]))
    print("TGT:", "".join([c if i in common_structure else '-' for i,c in enumerate(tgt_seq.sequence) ]))

    # Trim input PDB containing only common_structure residues
    pdb_pos2align_pos, align_pos2pdb_pos = translation_dicts(seq.sequence, modomics_from_pdb(pdb).sequence)

    print("MOD:", modomics_from_pdb(pdb))
    print(pdb_pos2align_pos)

    pdb_pos_trim = [align_pos2pdb_pos[c] for c in common_structure if c in align_pos2pdb_pos.keys()]
   
    print(pdb_pos_trim)
    print("MOD:", "".join([cc if ii + 1 in pdb_pos_trim else '-' for ii, cc in enumerate(modomics_from_pdb(pdb).sequence) ]))

    #exit()
    print("Eventual sequence length will be", len([c for c in tgt_seq.sequence if c != '-']))
    trim_pdb_to_res(pdb, pdb.replace('.pdb', '_trimmed_{}.pdb'.format(old_pdb)), pdb_pos_trim)

    # Thread them to their new identities.
    tgt_seq_for_thread = "".join([tgt_seq.sequence[c] for c in common_structure])

    command = [exe("rna_thread_and_minimize"),
        "-s", pdb.replace('.pdb', '_trimmed_{}.pdb'.format(old_pdb)), 
        "-seq", tgt_seq_for_thread.replace('-',''), 
        "-input_sequence_type", "MODOMICS", 
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
        "-guarantee_no_DNA", "true",
        "-ignore_zero_occupancy", "false",
        "-overwrite",
        "-set_weights", "other_pose", "0.0", "intermol", "0.0", "loop_close", "0.0", "free_suite", "0.0", "free_2HOprime", "0.0"]
    #print(command)
    subprocess.run(command)
    
    # Remove LINK records. They will be misleading because of a bug in the threading code.
    cull_LINKs(pdb.replace('.pdb', '_trimmed_{}_0001.pdb'.format(old_pdb)), "{}_input.pdb".format(old_pdb))

    old_seq2new_seq: Dict[int, int] = seq_map(seq.sequence, tgt_seq.sequence)

    # We would have to check that this is robust to any starting numbering (T:10-85)
    new_trimmed = [old_seq2new_seq[p] for p in pdb_pos_trim]

    numbering = dashed(new_trimmed, 'A')
    #print(numbering)
    
    # renumber input if needed? For sure to new chain.
    # ugly -- but numbering has spaces in it
    command = ["renumber_pdb_in_place.py", "{}_input.pdb".format(old_pdb)]
    command.extend(numbering.split())
    #print(command)
    subprocess.run(command)

    # Create the overall fasta.
    write_fasta("target.fasta", tgt_seq)

    # Do denovo run.
    command = [exe("rna_denovo"),
        "-s", "{}_input.pdb".format(old_pdb), 
        "-fasta", "target.fasta", 
        "-minimize_rna", "true",
        "-include_neighbor_base_stacks", "true",
        "-motif_mode", "true",
        "-score:weights", "stepwise/rna/rna_res_level_energy7beta.wts",
        "-set_weights", "other_pose", "0.0", "intermol", "0.0", "loop_close", "0.0", "free_suite", "0.0", "free_2HOprime", "0.0",
        "-guarantee_no_DNA", "true",
        "-use_legacy_job_distributor", "true",
        "-new_fold_tree_initializer", "true",
        "-nstruct", "1",
        "-out:file:silent", "{}_based_modeled.out".format(old_pdb)]
    
    if defer:
        with open("README_FARFAR_{}".format(old_pdb), "w") as f:
            f.write("{}\n\n".format(" ".join(command)))
    else:
        subprocess.run(command)
        subprocess.run(["extract_lowscore_decoys.py", "{}_based_modeled.out".format(old_pdb)])

