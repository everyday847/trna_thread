from thread_denovo_trna import dashed, simple_match, Sequence
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

def add_dash_recursive(template, trial, dashes, current_best):
    def filled_trial_seq(trial: str, ii: int) -> str:
        return trial[:ii]+'-'+trial[ii:]+'-'*(len(template)-len(trial)-1)

    #print("Call to add_dash_recursive with parameters:")
    #print(template.sequence)
    #print(trial)
    #print(dashes)
    n_left_to_add = len(template)-len(trial)
    print("Must add", n_left_to_add, "dashes")
    # Construct
    if len(dashes) == n_left_to_add:
        # We have enough.
        complete_trial_string = "-"*len(template)
        trial_index = 0
        for complete_index in range(len(template)):
            print(dashes)
            if complete_index in dashes: continue
            else:
                complete_trial_string = complete_trial_string[:complete_index]+trial[trial_index]+complete_trial_string[complete_index+1:]
                trial_index += 1
        print(complete_trial_string)
        score = simple_match(template, Sequence(complete_trial_string))
        print(current_best)
        if current_best is None or score > current_best[1]:
            current_best = (complete_trial_string, score)
        print(current_best)
    for ii in range(len(template)):
    	if ii in dashes: continue
    	else:
            new_dashes = list(dashes)
            new_dashes.append(ii)
    	
            current_best = add_dash_recursive(template, trial, new_dashes, current_best)

    #for ii in range(first_pos, len(template)-1):
    #    print(template.sequence, filled_trial_seq(trial, ii))
    #    score = simple_match(template, Sequence(filled_trial_seq(trial, ii)))
    #    if current_best is None or score > current_best[1]:
    #        current_best = (filled_trial_seq(trial, ii), score)
    #    print(template.sequence, filled_trial_seq(trial, ii), score)
    #    if len(trial[:ii]+'-'+trial[ii:]) < len(template) \
    #        and trial[ii:] != '-'*(len(trial)-ii-1): 
    #        current_best = add_dash_recursive(template, filled_trial_seq(trial, ii), ii, current_best)
    return current_best

def test_add_dash_recursive():
    print("\n")
    #add_dash_recursive(Sequence('-A--A-'), 'AA', 0)
    current_best = ('', -100000)
    best_sequence = add_dash_recursive(Sequence('-A--A-'), 'AA', [], current_best)

