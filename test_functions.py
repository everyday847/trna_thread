from lib.str_manip import dashed
from thread_denovo_trna import simple_match
from lib.Sequence import Sequence
import mypy
from typing import List, Dict, Tuple

def test_dashed():
    """
    We're okay with a terminating space because we can just strip it off if it
    ever matters.
    """
    assert(dashed([1,2,3,7,8,9], 'A') == "A:1-3 A:7-9 ")

def test_dashed_solo():
    """
    We're okay with a terminating space because we can just strip it off if it
    ever matters.
    """
    assert(dashed([1,2,3,7,8,9,15], 'A') == "A:1-3 A:7-9 A:15 ")

def test_dashed_solo_mid():
    """
    We're okay with a terminating space because we can just strip it off if it
    ever matters.
    """
    assert(dashed([1,2,3,5,7,8,9], 'A') == "A:1-3 A:5 A:7-9 ")

def test_simple_match_identity():
    """
    aaa
    """
    assert(simple_match(Sequence('A-'), Sequence('A-')) == 10)
    assert(simple_match(Sequence('G-'), Sequence('G-')) == 10)
    assert(simple_match(Sequence('GGAA---'), Sequence('GGAA---')) == 35)

def test_simple_match_gap_penalty():
    assert(simple_match(Sequence('A-'), Sequence('-A')) == -20)

def test_simple_match_purine_eq():
    assert(simple_match(Sequence('A-'), Sequence('G-')) == 7)

def test_simple_match_pyrimidine_eq():
    assert(simple_match(Sequence('U-'), Sequence('C-')) == 7)

def test_simple_match_mod_U_eq():
    assert(simple_match(Sequence('U-'), Sequence('T-')) == 8)
