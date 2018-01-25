import pickle
import re

extractor = re.compile("(chr[0-9|XY]{1,2}):([\d]+)-([\d]+)\_\[\'(-?1)")

def parse_kv(key, value, extractor):
     guide = key[:-3]
     string_val = value[0]
     match = extractor.match(string_val)
     if match:
             chrom, start, end, strand = match.groups()
             return guide, chrom, start, end, strand
     else:
             return guide, string_val, "error", "error", "error"

with open('gRNA_target_coordinates_annotation_dict.pkl','r') as f:
	d = pickle.load(f)

with open('guide_to_target_region.csv','w') as f:
	for k in d.keys():
		guide, chrom, start, end, strand = parse_kv(k, d[k], extractor)
		print >> f, ','.join([guide, chrom, start, end, strand])
