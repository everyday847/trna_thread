import os, mypy
from typing import Dict, List, Tuple
from lib.Sequence import Sequence
import numpy as np

###
#   Sequence translation (among Rosetta and Modomics nomenclatures)
###

# yW won't work for some reason? change to YYG always. ugh.
# YG is yW, AET is MTA, G7M is 7MG, QUO is   Q, 1RN is CNS, YYG is yW [ish], 70U is MST, 12A is M26, 2MU is 52U, 6IA is I6A
# oh crap: the first stuff breaks mod_to_tlc for situations where we disagree with PDB. This is a whole mess for
# YYG.
tlc_to_mod = { #                                                                       note this is the enantiomer at the AA
               "GTP": 'G', "YYG": 'Y', "AET": 'E', "G7M": '7', "QUO": 'Q', "1RN": '$', "YYG": 'Y', "  I": 'I',
               "70U": '3', "12A": '[', "2MU": '\\', "6IA": '+',

               "  A": 'A', "  C": 'C', "  G": 'G', "  U": 'U', "1MA": '"', "2MA": '/', "I6A": '+', "MIA": '*',
               "MA6": '=', "T6A": '6', "AET": 'E', "M26": '[', "OMA": ':', "INO": 'I', "M1I": 'O', "RIA": '^',
               "IHA": '`', "S2C": '%', "OMC": 'B', "A4C": 'M', "5MC": '?', "3MC": '\'', "K2C": '}', "F5C": '>',
               "52C": '<', "1MG": 'K', "2MG": 'L', "OMG": '#', "M2G": 'R', "M3G": '|', "7MG": '7', "FA7": '(',
               "  Q": 'Q', "YYG": 'Y', "O2W": 'W', "2SU": '2', "OMU": 'J', "4SU": '4', "5MU": 'T', "2ST": 'F',
               "52U": '\\', "NMT": '{', "NST": 'S', "CST": '&', "CMT": '~', "MOT": '1', "MST": '3', "OAU": 'V',
               "MOU": '5', "CNT": '!', "CNS": '$', "CNM": ')', "APU": 'X', "MHU": ',', "H2U": 'D', "PSU": 'P',
               "1PU": ']', "MPU": 'Z' }
mod_to_tlc: Dict[str, str] = dict(zip(tlc_to_mod.values(), tlc_to_mod.keys()))

U_equivs = {'T', 'P', '$', '3', '\\', '2', 'J', '4', 'F', '{', 'S', '&', '~', '1', 'V', '5', '!', '$', ')', 'X', ',', 'D', ']', 'Z'}
A_equivs = {'[', '+', '"', '/', '*', '=', '6', 'E', ':', '^', '`'}
C_equivs = {'%', 'B', 'M', '?', '\'', '}', '>', '<'}
G_equivs = {'7', 'Y', 'Q', 'I', 'O', 'K', 'L', '#', 'R', '|', '(', 'W'}

def mod_from_tlc(tlc: str) -> str:
    """
    Translates modomics single letter codes from 3LCs
    """

    if tlc in tlc_to_mod:
        return tlc_to_mod[tlc]

    print("Unrecognized", tlc)
    exit()

def modomics_from_pdb(pdb: str) -> Sequence:
    """
    Doesn't pair with secstruct (yet). dssr?
    """

    pdblines = [] # type: List[str]
    modomics_seq = ""
    with open(pdb) as f:
        pdblines = f.readlines()
    for l in pdblines:
        if " C4'" in l: modomics_seq += mod_from_tlc(l[17:20])

    # Replace with Rosetta-native SS determination?

    return Sequence(modomics_seq, '.'*len(modomics_seq))

def ann_to_mod(ann_seq: str) -> str:
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

def annotated_seq_of(mod_seq: str) -> str:
    ann_seq = ""
    for c in mod_seq:
        if c == 'A': ann_seq += 'a'
        elif c == 'C': ann_seq += 'c'
        elif c == 'G': ann_seq += 'g'
        elif c == 'U': ann_seq += 'u'
        else: ann_seq += 'X[' + mod_to_tlc[c] + ']'
    return ann_seq

def dashed(positions: List[int], chain: str) -> str:
    """
    Take a set of positions (i.e., 1 2 3 7 8 9) and a chain ("A"), and produce a
    repr like "A:1-3 A:7-9"
    """
 
    #print(positions)
    positions = sorted(positions)
    acc_start = None
    acc_end = None
    numbering = ""
    for p in range(1, max(positions)+1):
        #print("examining", p)
        if p in positions and p != max(positions):
            #print("\tfound in positions", p)
            if acc_start is None:
                acc_start = p
            acc_end = p
            #print("now acc_start acc_end", acc_start, acc_end)
        elif p != max(positions):
            if acc_start is not None and acc_end is not None:
                if acc_start != acc_end:
                    numbering += "{ch}:{start}-{end} ".format(ch=chain, start=acc_start, end=acc_end)
                else:
                    numbering += "{ch}:{only} ".format(ch=chain, only=acc_start)
                acc_start = None
                acc_end = None
        else:
            if p in positions:
                acc_end = p
            if acc_start is None:
                numbering += "{ch}:{only} ".format(ch=chain, only=acc_end)
            else:
                numbering += "{ch}:{start}-{end} ".format(ch=chain, start=acc_start, end=acc_end)

    return numbering

def all_bp_partners(ss_str: str) -> Dict[int, int]:
    """
    We push back on a dict-list (per char type), then when we find a
    complement we pop -- giving a pair of popped + current.
    """
    final = {} # type: Dict[int, int]
    paired = {'(': [], '[': [], '{': []} # type: Dict[str, List[int]]
    compl = {'[': ']', ']': '[', '(': ')', ')': '(', '{': '}', '}': '{'}
    for i, c in enumerate(ss_str):
        if c == '.': continue
        if c in paired.keys(): paired[c].append(i)
        if c == ')' or c == ']' or c == '}':
            pair = paired[compl[c]].pop()
            final[pair] = i
            final[i] = pair
    return final
