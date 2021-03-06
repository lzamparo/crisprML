__author__ = 'Alexendar Perez'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""compute specificity score, Hamming, and Levinstein distance neighborhoods for strings"""

#################
#               #
#   Libraries   #
#               #
#################

import sys
import os
import pickle
import argparse
import sqlite3
import gzip

import numpy as np
import pandas as pd
from Bio import trie


#############################
#							#
#	CFD Scoring Functions	#
#							#
#############################

def calc_cfd(wt,sg,pam,mm_scores,pam_scores):
	#mm_scores,pam_scores = get_mm_pam_scores()
	score = 1
	sg = sg.replace('T','U')
	wt = wt.replace('T','U')
	s_list = list(sg)
	wt_list = list(wt)
	for i,sl in enumerate(s_list):
		if wt_list[i] == sl:
			score*=1
		else:
			try:
				key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
				score*= mm_scores[key]
			except KeyError:
				continue
	score*=pam_scores[pam]
	return (score)

def get_mm_pam_scores(mms,pams):
	try:
		mm_scores = pickle.load(open(mms,'rb'))
		pam_scores = pickle.load(open(pams,'rb'))
		sys.stdout.write('CFD scoring matrices loaded\n')
		return (mm_scores,pam_scores)
	except:
		raise Exception("Could not find file with mismatch scores or PAM scores")

def revcom(s):
	basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
	letters = list(s[::-1])
	letters = [basecomp[base] for base in letters]
	return ''.join(letters)

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--infile',help='absolute filepath to file with gRNA sequences',required=True)
	parser.add_argument('-o','--outdir',help='absolute filepath to output directory',required=True)
	parser.add_argument('-k','--kmer',help='absolute filepath to kmers_counted.txt file',required=True)
	parser.add_argument('-t','--trie',help='absolute filepath to trie.dat file',required=True)
	parser.add_argument('-m','--mismatch',help='absolute filepath to mismatch_score.pkl for CFD',required=True)
	parser.add_argument('-p','--pam',help='absolute filepath to pam_scores.pkl for CFD',required=True)
	parser.add_argument('--header',help='boolian value of whether header is present in infile, default = True',default=True)
	parser.add_argument('--sequence_field',help='if sequences not in first field of file, default = 0',default=0)
	parser.add_argument('--cpf1',help='cpf1 enzyme processing',default=False)

	args = parser.parse_args()
	in_file = args.infile
	outdir = args.outdir
	kmer_counts_file = args.kmer
	trie_file = args.trie
	mismatch_score = args.mismatch
	pam_score = args.pam
	header = args.header
	sequence_field = args.sequence_field
	cpf1 = args.cpf1

	return in_file,outdir,kmer_counts_file,trie_file,mismatch_score,pam_score,header,int(sequence_field),cpf1

def load_pickle(infile):
	"""load pickle file

	:param infile: absolute filepath to pickle
	:return: deserialized pickle

	"""
	with open(infile, 'r') as in_file:
		data = pickle.load(in_file)

	return data

def sequence_data_extraction(data,header,sequence_field):
	""" extract sequence data from data array

	:param data: numpy array of data, first output of sequence_file_read_in()
	:param header: boolian value indicating if header is present in data array; will skip first line if header present
	:param sequence_field: if sequence data is not in 0 field specify; defaults to 0 field
	:return: sequence data column

	"""
	if header:
		sys.stdout.write('skipping first line in data array due to header\n')
		if sequence_field:
			sys.stdout.write('sequence field specified as %s\n' % sequence_field)
			sequence_data = data[1:,sequence_field]
		else:
			sys.stdout.write('sequence field defaulted to 0\n')
			sequence_data = data[1:,0]
	else:
		if sequence_field:
			sys.stdout.write('sequence field specified as %s\n' % sequence_field)
			sequence_data = data[:,sequence_field]
		else:
			sys.stdout.write('sequence field defaulted to 0\n')
			sequence_data = data[:,0]

	return sequence_data

def sequence_file_read_in(in_file):
	"""read in file with sequences like gRNAs

	:param in_file: absolute filepath to file containing sequences
	:return: numpy array representation of data accessed through either pickle or pandas modules

	"""
	sys.stdout.write(
		'%s is being used to compute features for ClassTask\n***Sequence data should be in first field***\n' % in_file)
	try:
		sys.stdout.write('attempting to read %s as pickle\n' % in_file)
		file_format = 'pickle'
		data = load_pickle(in_file)
	except:
		try:
			sys.stdout.write('attempting to read %s with pandas as excel\n' % in_file)
			file_format = 'excel'
			data = np.array(pd.read_excel(in_file, header=None))
		except:
			try:
				sys.stdout.write('attempting to read %s with pandas as text file\n' % in_file)
				file_format = 'text'
				data = np.array(pd.read_table(in_file, header=None))
			except:
				sys.stderr.write('%s file format not recognized as pickle, excel, or text file; aborting\n' % in_file)
				sys.exit(1)

	sys.stdout.write('%s successfully read in as %s\n' % (in_file, file_format))
	return data

def load_trie(trie_file):
	"""deserialize trie

	:param trie_file: serialized trie from BioPython and produced through GuideScan processer.py: x__all_trie.dat
	:return: deserialized trie object

	"""
	tr_file = open(trie_file, 'r')
	tr = trie.load(tr_file)
	tr_file.close()
	sys.stdout.write('trie loaded into memory from %s \n' % trie_file)

	return tr

# TODO: this is better stored as a DB; it also does not get used as the docstring 
# suggests: they keys are parts[1], which are *genomic positions*, and the values
# are kmers.
def kmer_exact_occurrence_dictionary(kmer_counts_file):
	"""generate genome-wide kmer occurrence dictionary as an sqlite database.  If the 
	filename suggests that the kmer_counts_file is already a .db, then just return the cursor.

	:param kmer_counts_file: absolute filepath to XXX_all_kmers_counted.txt.gz file, or a .db file
	:return: connection object to the kmer-count dictionary 

	"""


	if kmer_counts_file.endswith('db'):
		conn = sqlite3.connect(kmer_counts_file)
		return conn
	

	kmer_dictionary = {}
	records = 0

	sqlite_file = 'hg38_kmers.db'
	table_name = 'kmer_counts'  
	first = 'kmer' 
	first_type = 'text'  
	second = 'count'
	second_type = 'INTEGER'

	# Connect to the database file
	try:
		conn = sqlite3.connect(os.path.join(os.path.dirname(kmer_counts_file),sqlite_file))
	except: 
		sys.stdout.write('Cannot open the sqlite db! \n')
		sys.exit()
		
	c = conn.cursor()

	# Create the table
	c.execute('CREATE TABLE {tn} ({fc} {ft}, {sc} {st})'\
        .format(tn=table_name, fc=first, ft=first_type, sc=second, st=second_type))

	my_open = gzip.open if kmer_counts_file.endswith('.gz') else open

	with my_open(kmer_counts_file, 'r') as kmers:
		for line in kmers:
			clean_line = line.lstrip().rstrip()
			parts = clean_line.split()
			if kmer_dictionary.has_key(parts[1]):
				sys.stdout.write('kmer duplication detected %s %s \n' % (parts[1], parts[0]))
			else:
				kmer_dictionary[parts[1]] = parts[0]

			records += 1

			# dump dict into the database, then reset it
			if records > 100000:
				for k, v in kmer_dictionary.items():
					c.execute("INSERT INTO kmer_counts VALUES (?,?)", (k, v))
				kmer_dictionary = {}
				records = 0

		# handle the remaining records
		for k, v in kmer_dictionary.items():
			c.execute("INSERT INTO kmer_counts VALUES (?,?)", (k, v))

		sys.stdout.write('kmer dictionary generated \n')
	
	# commit the changes to the db
	conn.commit()
	return conn

def add_features_to_feature_array(feature_array,augmenting_array):
	"""add new features to features array

	:param feature_array: numpy array with sequences as UI and previously computed features
	:param augmenting_array: numpy array with new features to be added to feature array
	:return: new feature array with features from the augmented array added

	"""
	equivilence = np.all(feature_array[:, 0].reshape(feature_array.shape[0], 1) == augmenting_array[:, 0].reshape(augmenting_array.shape[0], 1))
	if equivilence:
		try:
			feature_array = np.concatenate((feature_array,augmenting_array[:,1:]),1)
			return feature_array
		except IndexError:
			feature_array = np.concatenate((feature_array.reshape(feature_array.shape[0],1), augmenting_array[:, 1:]), 1)
			return feature_array

	else:
		sys.stdout.write('original data array and new features array NOT in same order: attempt sort\n')
		feature_array = feature_array[feature_array[:, 0].argsort()]
		augmenting_array = augmenting_array[augmenting_array[:, 0].argsort()]
		return add_features_to_feature_array(feature_array,augmenting_array)

def hamming_distance(s1, s2):
	"""calculate the Hamming distance between two bit strings

	:param s1: first string
	:param s2: second string
	:return: Hamming distance between s1 and s2

	"""
	assert len(s1) == len(s2)
	return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def query_db(c, key):
	"""
	Query the kmer file sqlite3 database with cursor `c` to extract the count associated with kmer `key`
	"""
	try:
		result = c.execute("SELECT count FROM kmer_counts WHERE kmer == ?", (key,))
		result_list = result.fetchall()
		if len(result_list) == 1:
			return result_list[0][0]
		elif len(result_list) == 0:
			return 0
	except:
		sys.stderr.write('querying db returned loads of hits for {0}: {1}'.format(key, ' '.join([str(r) for r in result_list])))
		

def compute_specificity_score_and_mismatch_neighborhoods(sequence_data, final_header, kmer_dictionary_cursor, tr, mm_scores,
														 pam_scores,cpf1):
	"""compute GuideScan based features

	:param sequence_data: numpy array with sequence data in 0 field
	:param final_header: numpy array with header information
	:param kmer_dictionary_cursor: cursor for the sqlite database; output of kmer_exact_occurrence_dictionary()
	:param tr: trie datastructure from load_trie() function
	:param mm_scores: first output of get_mm_pam_scores()
	:param pam_scores: second output of get_mm_pam_scores()
	:return: feature array with GuideScan derived features

	"""
	distance = 3
	guidescan_array, seq_array = np.zeros((sequence_data.shape[0], 10)), np.empty((sequence_data.shape[0], 1)).astype(
		str)
	for j, on_target_sequence in enumerate(sequence_data[:, 0]):

		# sequence array value
		on_target_sequence_value = on_target_sequence

		print "Processing guide " + str(j)
		sys.stdout.flush()

		# guidescan format
		if cpf1:
			on_target_sequence = 'TTTN%s' % (on_target_sequence)
		else:
			on_target_sequence = '%sNGG' % (on_target_sequence)

		# query trie, get all near matches
		query_sequences = tr.get_approximate(on_target_sequence, distance)

		# specificity score lists
		cfd_lst, writeout_lst = [], []

		# neighborhood enumeration dictionaries
		hamming_key, hamming_distance_dict, levinstein_key, levinstein_distance_dict = {}, {0: 0, 1: 0, 2: 0,
																							3: 0}, {}, {1: 0, 2: 0,
																										3: 0}

		for i in query_sequences:
			# occurrence of sequence in genome
			ot_sequence_occurence = int(query_db(kmer_dictionary_cursor,i[0]))

			if hamming_distance(on_target_sequence, i[0]) <= distance:
				# record key
				if hamming_key.has_key(i[0]):
					continue
				else:
					if cpf1:
						pass
					else:
						# cfd computation
						pam = i[0][-2:]
						sg = i[0][:-3]
						cfd_score = calc_cfd(on_target_sequence, sg, pam, mm_scores, pam_scores)
						total_cfd_contribution = cfd_score * float(ot_sequence_occurence)
						cfd_lst.append(total_cfd_contribution)

					# augment count for Hamming neighbors at n mismatches
					if hamming_distance_dict.has_key(hamming_distance(on_target_sequence, i[0])):
						hamming_distance_dict[hamming_distance(on_target_sequence, i[0])] = \
							hamming_distance_dict[hamming_distance(on_target_sequence, i[0])] + ot_sequence_occurence
						hamming_key[i[0]] = 1

					# establish count for Hamming neighbors at n mismatches
					else:
						hamming_distance_dict[hamming_distance(on_target_sequence, i[0])] = ot_sequence_occurence
						hamming_key[i[0]] = 1

			else:
				# record key
				if levinstein_key.has_key(i[0]):
					continue
				else:
					# augment count for Levinstein neighbors at n mismatches
					if levinstein_distance_dict.has_key(i[2]):
						levinstein_distance_dict[i[2]] = levinstein_distance_dict[i[2]] + ot_sequence_occurence
						levinstein_key[i[0]] = 1

					# establish count for Hamming neighbors at n mismatches
					else:
						levinstein_distance_dict[i[2]] = ot_sequence_occurence
						levinstein_key[i[0]] = 1

		# cfd composite specificity score
		if cpf1:
			cfd_aggregate_score = 0
		else:
			cfd_array = np.array(cfd_lst)
			cfd_aggregate_score = 1.0 / (cfd_array.sum())

		# fill in features into feature array
		seq_array[j, 0] = on_target_sequence_value
		guidescan_array[j, 0] = cfd_aggregate_score
		guidescan_array[j, 1] = int(hamming_distance_dict[0])
		guidescan_array[j, 2] = int(hamming_distance_dict[1])
		guidescan_array[j, 3] = int(hamming_distance_dict[2])
		guidescan_array[j, 4] = int(hamming_distance_dict[3])
		guidescan_array[j, 5] = sum(hamming_distance_dict.values())
		guidescan_array[j, 6] = int(levinstein_distance_dict[1])
		guidescan_array[j, 7] = int(levinstein_distance_dict[2])
		guidescan_array[j, 8] = int(levinstein_distance_dict[3])
		guidescan_array[j, 9] = sum(levinstein_distance_dict.values())

		"""
		sys.stdout.write('Hamming enumerated neighbors = %s\nLevinstein enumerated neighbors = %s\nSpecificity score = %s\nhamming sequence neigbors = %s\n'
						 % (sum(hamming_distance_dict.values()),sum(levinstein_distance_dict.values()),cfd_aggregate_score,len(cfd_lst)))
		"""

	# generate final augmented features array
	seq_guidescan_array = np.concatenate((seq_array, guidescan_array), 1)

	sequence_data = add_features_to_feature_array(sequence_data, seq_guidescan_array)
	header_value = np.array(['Specificity_Score', 'Occurrences_at_Hamming_0', 'Occurrences_at_Hamming_1',
							 'Occurrences_at_Hamming_2', 'Occurrences_at_Hamming_3', 'Sum_Hamming_Neighbors',
							 'Occurrences_at_Levinstein_1', 'Occurrences_at_Levinstein_2',
							 'Occurrences_at_Levinstein_3',
							 'Sum_Levinstein_Neighbors']).reshape(1, 10)
	final_header = np.concatenate((final_header, header_value), 1)
	sys.stdout.write('GuideScan based features computed\n')
	return sequence_data, final_header

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():

	"""
	in_file = '/Users/pereza1/Projects/Jo/data/gecko_proper_excel/mouse_library_A_gecko.xlsx'
	header = True
	sequence_field = 0
	"""

	# user inputs
	in_file,outdir,kmer_counts_file,trie_file,mismatch_score,pam_score,header,sequence_field,cpf1 = arg_parser()


	# data read in
	data = sequence_file_read_in(in_file)

	# sequence data extraction
	sequence_data = sequence_data_extraction(data,header,sequence_field)
	sequence_data = sequence_data.reshape(sequence_data.shape[0],1)
	final_header = np.array(['sequence']).reshape(1,1)

	# compute or load kmer dictionary object
	kmer_dictionary = kmer_exact_occurrence_dictionary(kmer_counts_file)
	kmer_dictionary_cursor = kmer_dictionary.cursor()

	# load CFD scoring matrices
	mm_scores, pam_scores = get_mm_pam_scores(mismatch_score, pam_score)

	# load trie
	tr = load_trie(trie_file)

	# compute specificity score and mismatch neighborhoods
	sequence_data,final_header = compute_specificity_score_and_mismatch_neighborhoods(sequence_data,final_header,
																					  kmer_dictionary_cursor,tr,mm_scores,
																					  pam_scores,cpf1)

	# generate final feature arrays
	final_feature_array = np.concatenate((final_header,sequence_data),0)
	#final_feature_array_standardized = np.concatenate((final_header,sequence_data_standardized),0)
	sys.stdout.write('final feature arrays generated\n')

	# write output to csv
	column_length = final_feature_array.shape[1]
	np.savetxt('%s/raw_features_computed_%s.csv' % (outdir,in_file.split('/')[-1].split('.')[0]), final_feature_array,
			   fmt='%' + '%ss' % (column_length), delimiter=',')

	#np.savetxt('%s/standarized_features_computed_%s.csv' % (outdir,in_file.split('/')[-1].split('.')[0]), final_feature_array_standardized,
	#		   fmt='%' + '%ss' % (column_length), delimiter=',')

	sys.stdout.write('final arrays written to csv\n%s\n' % ('%s/features_computed_%s.csv' % (outdir,in_file.split('/')[-1].split('.')[0])))

	# close the kmer_dictionary db
	kmer_dictionary.close()
	
	# completion stdout
	sys.stdout.write('feature generation for %s complete\n' % (in_file))

if __name__ == '__main__':
	main()