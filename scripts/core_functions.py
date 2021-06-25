def revcom(s):
	"""
	:param s: str
	:return: reverse complemented string
	"""
	basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V',
                'V':'B','H':'D','D':'H','-':'-'}
	letters = list(s[::-1])
	letters = [basecomp[base] for base in letters]
	return ''.join(letters)
