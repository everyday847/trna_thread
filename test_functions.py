from thread_denovo_trna import dashed

def test_dashed():
	"""
	We're okay with a terminating space because we can just strip it off if it
	ever matters.
	"""
	assert(dashed([1,2,3,7,8,9], 'A') == "A:1-3 A:7-9 ")