from thread_denovo_trna import dashed, simple_match

def test_dashed():
    """
    We're okay with a terminating space because we can just strip it off if it
    ever matters.
    """
    assert(dashed([1,2,3,7,8,9], 'A') == "A:1-3 A:7-9 ")

def test_simple_match_identity():
    """
    aaa
    """
    assert(simple_match('A-', 'A-') == 10)
    assert(simple_match('G-', 'G-') == 10)
    assert(simple_match('GGAA---', 'GGAA---') == 35)

def test_simple_match_gap_penalty():
    assert(simple_match('A-', '-A') == -20)

def test_simple_match_purine_eq():
    assert(simple_match('A-', 'G-') == 7)

def test_simple_match_pyrimidine_eq():
    assert(simple_match('U-', 'C-') == 7)

def test_simple_match_mod_U_eq():
    assert(simple_match('U-', 'T-') == 8)