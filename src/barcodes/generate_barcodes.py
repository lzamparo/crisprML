import sys

def allstrings(alphabet, length, limit):
    """Find the list of all strings of 'alphabet' of length 'length'"""
    
    if length == 0: return []
    
    c = [[a] for a in alphabet[:]]
    if length == 1: return c
    
    c = [[x,y] for x in alphabet for y in alphabet]
    if length == 2: return c
    
    for l in range(2, length):
        c = [[x]+y for x in alphabet for y in c]
	if len(c) > limit:
		break        

    return c

if __name__ == "__main__":
    
	length = sys.argv[1]
	alphabet = ['A','C','G','T']
	limit = sys.argv[2]
	outfile = sys.argv[3]

	my_strings = allstrings(alphabed, length, limit)
	with open(outfile,'w') as barcodes:
		for s in my_strings:
			print(s, file=barcodes)	
