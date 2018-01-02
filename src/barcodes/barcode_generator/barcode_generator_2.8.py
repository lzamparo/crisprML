#!/usr/local/bin/python

'''BARCODE GENERATOR
    BY Luca Comai and Tyson Howell
    December 2009'''

print '\n\nBARCODE GENERATOR by LC and TH\n\t--***---\nPlant Biology and Genome Center\n\tUC Davis\n'
print '''This program generates barcodes
of a desired length, distance, and GC content
The primer sequences currently used
are the illumina paired end primers\n'''


#___________________________________________________raw input

print 'Enter LENGTH as an integer (i.e. 4)'

# ask the desired number of bases in barcode
length = int(raw_input("Barcode length: "))

print '\nEnter the number of barcodes (default is LENGTH x 5)'

# ask how many barcodes should be made
number = raw_input("Total number of barcodes: ")
if number == '':
    number = length*5
else:
    number = int(number)

print '\nEnter the minimum number of different bases between barcodes (default is LENGTH/2, i.e. 7->3, 4->2)'

# ask what is the least number of bases that must be different
# between any two barcodes
diffs = raw_input("Min. no. of different bases between barcodes: ")
if diffs == '':
    diffs = length//2
else:
    diffs = int(diffs)

print '\nEnter desired GC content range in percentages (i.e. 50 ->50%)'

# ask for desired GC content range
mingc = raw_input("Minimum GC content (default is 0):")
if mingc == '':
    mingc = 0
else:
    mingc = float(mingc) / 100

maxgc = raw_input("Maximum GC content (default is 100):")
if maxgc == '':
    maxgc = 100
else:
    maxgc = float(maxgc) / 100

print '\nThe default number of random codes to test is 10000. Do not enter more than a million'

# ask what is the maximum number of random codes to be tested
attempts = raw_input('How many attempts?:')
if attempts == '':
    attempts = 10000
else:
    attempts = int(attempts)

#___________________________________________________process

# make list of the four bases
l1 = ['a', 'c', 'g', 't']

# initialze the barcode list
barcode_list = []

# initialize the first barcode
first_barcode = []

# prime the tested list, for future counting
tested = []

# function to determine GC content
def gc_cont (bar_code):
    gc = 0.0
    for base in range(length):
        if bar_code[base] == 'c' or bar_code[base] == 'g':
            gc += 1
        else:
            gc += 0
    cont = gc / length
    return cont

# import random module
import random

# make the first barcode
# add first barcode to barcode list. This is needed for the
# first comparison of "compare_barcode" function
while barcode_list == []:
    for i in range(length):
        first_barcode.append(random.choice(l1))
    if gc_cont(first_barcode) <= maxgc and gc_cont(first_barcode) >= mingc:
        barcode_list.append(first_barcode)
    else:
        first_barcode = []
        

# the barcode "cradle": a place where each barcode will sit
barcode = []

#___________________________________________________define functions

# function makes the barcode
def make_barcode(length):
    global barcode
    # empties the barcode cradle
    barcode = []
    for i in range(length):
        barcode.append(random.choice(l1))

# barcode is tested vs the previously generated barcodes
def compare_barcode(length, barcode_l):
    count = 0
    global barcode
    # run barcode creator
    make_barcode(length)
    # keep track of it
    tested.append(barcode)
    # testing of barcode
    if barcode not in barcode_list:
        global count_list
        count_list = []
        # compare to barcodes in list
        for bc in barcode_l:
            # matches to existing barcodes
            # are scored as points
            count = 0
            for pos in range(length):
                if barcode[pos] == bc[pos]:
                    count += 1
                else:
                    count += 0
            # for each barcode a list of scores is made
            count_list.append(count)
        # if the barcode has enough unique bases
        # and the proper GC content, it is added
        # to the list of good barcodes
        if max(count_list) > length-diffs:
            count_list = []
        elif gc_cont(barcode) <= maxgc and gc_cont(barcode) >= mingc:
            barcode_list.append(barcode)
            count_list = []
        else:
            count_list = []
    else:
        pass
               
#___________________________________________________run functions
    
# initialize count

count_list = []

# program stalls if too many attempts are allowed
# and few barcodes remain to be discovered
# this loop keeps the attempts within the range allowed
while len(tested) < attempts:
    if len(barcode_list) < number:
        compare_barcode(length, barcode_list)
    else:
        break

barcode_list.sort()


print "\n\nRESULTS\n\ngood barcodes and GC content:"

for i in barcode_list:
    print i, int(gc_cont(i)*100), '%'

print '\nnumber of tested barcodes:'

print len(tested)

print '\nnumber of good barcodes:'

print len(barcode_list)

#count base composition in each of the barcode position
from collections import defaultdict

print '\nbase compositions by position'

for pos in range(length):
    list_l = [i[pos] for i in barcode_list]
    base_count = defaultdict(int)
    for base in list_l:
        base_count[base] += 1
    print base_count.keys(), base_count.values()

#___________________________________________________manipulate barcodes and print results

# make a file for the barcoded primer sequence
# open file 
barfile = open('barcode.txt', 'wb')

print >> barfile, 'These are the barcodes of length '+str(length)+' with a distance of '+str(diffs)+' bases\n'

#___________________________________________________primer sequence

# modify the primer as desired. This seq is the
# illumina PE primer
primerA = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'

# note that the adapter is originally as below
# GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
# however, an A is added to allow the sequencing
# primer to anneal

# this is the complementary primer from illumina
# same for regular or PE 
primerB = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'

# the output will be:
# >adA2_cccca
# ccccaAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
# >adB2_cccca
# ACACTCTTTCCCTACACGACGCTCTTCCGATCTtggggT
# note that barcode is in lower case

# define adapter names
name_rootA = '>adA2_'

name_rootB = '>adB2_'

# initialize a holder name for the complement of barcode
comp_barcode = ''

# define function to derive complement of any seq
# try/except/finally is to make sure that the 'maketrans'
# has been imported
def reverse_comp(seq):
    try:
        maketrans  
    except NameError:
        from string import maketrans
    finally:
        comp_table = maketrans('actg','tgac')
        global comp_barcode
        comp_barcode = seq[::-1].translate(comp_table)

for i in barcode_list:
    j = ''.join(i)
    reverse_comp(j)
    print >> barfile, '%s%s\n%s%s\n%s%s\n%s%sT\n\n' %  (name_rootA, j, j, primerA, name_rootB, j, primerB, comp_barcode)

print '''\nA file called barcode.txt has been generated.
It contains the adapter sequences with each
barcode in lower case\n\n'''

# close file or it does not get updated
barfile.close()

#___________________________________________________log

# changes from 2.7:
# rephrase raw input queries
# default barcodes to test to 10,000
# clean up shell results

