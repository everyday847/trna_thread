import mypy, os
from typing import List, Tuple
from lib.Sequence import Sequence

###
#   Data file manipulation: reading in sequence data and PDB file references.
###

def exe(name: str) -> str:
    """
    For distribution, once I've merged into master: look for ROSETTA in 
    os.environ and return exe path based on that.
    """
    import os
    if os.path.exists("/home/andy"):
        return "/home/andy/Rosetta/main/source/bin/" + name
    else:
        return "/Users/amw579/dev_ui/Rosetta/main/source/bin/" + name

def my_loc() -> str:
    """
    Since moving this to a library file... now it gives an undesired path.
    """
    import os
    return os.path.dirname(os.path.realpath(__file__)).replace('lib', '')

def get_seqs_mfa(fn: str) -> List[Tuple[str, Sequence]]:
    """
    Get sequences out of the tRNAdb mfa
    """
    lines = [] # type: List[str]
    with open(fn) as f:
        lines = f.readlines()
    sequences = [('', Sequence(l.strip())) for l in lines[1::2]]
    return sequences

def get_seqs(fn: str) -> List[Tuple[str, Sequence]]:
    lines = [] # type: List[str]
    with open(fn) as f:
        lines = f.readlines()
    #sequences = {l.strip().split()[0]: l.strip().split()[1] for l in lines if "ALIGNED_TO" in l }
    sequences = [(l.strip().split()[0], Sequence(l.strip().split()[1]))
        for l in lines if "PDB_SEQ" in l]
    return sequences
