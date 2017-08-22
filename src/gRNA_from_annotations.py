__author__ = 'Alexendar Perez'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""select gRNAs from a set that meet certain annotation requirements"""

#################
#               #
#   Libraries   #
#               #
#################

import sys
import pickle
import argparse
from collections import defaultdict

from bx.intervals.intersection import IntervalTree

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--annotations',help='absolute filepath to annotation BED file',required=True)
    parser.add_argument('-i','--infile',help='absolute filepath to GuideScan *_all_kmers.txt file',required=True)
    parser.add_argument('-o','--outdir',help='absolute filepath to output directory',required=True)
    parser.add_argument('--cpf1',help='if Cpf1 infile then --cpf1 = 1, else process as Cas9, default = 0',default=0)

    args = parser.parse_args()
    annotations = args.annotations
    gRNA_file = args.infile
    outdir = args.outdir
    cpf1 = args.cpf1

    return annotations,gRNA_file,outdir,int(cpf1)

def create_interval_tree(annotations):
    """create annotation tree datastructure

    :param annotations: absolute filepath to annotation BED file
    :return: interval tree

    """
    with open(annotations,'r') as infile:
        interval_tree = create_inttree_from_file(infile)
        sys.stdout.write('interval tree for %s generated\n' % annotations)

    return interval_tree

def gRNA_w_annotations_writeout(target_sequence_overlapping_exon, outdir, cpf1):
    """write out to file those gRNAs which fulfill annotation requirement

    :param target_sequence_overlapping_exon: list object, first output object of gRNA_annotation_through_interval_tree()
    :param outdir: absolute filepath to output directory
    :param cpf1: integer value, 0 indicates Cpf1 enzyme is being passed, != 0 means Cas9 is being passed
    :return: output file that can be passed to specificity_score_distance_neighbors()

    """
    if cpf1 == 1:
        enzyme = 'Cpf1'
    else:
        enzyme = 'Cas9'

    unique_dictionary = {}
    with open('%s/%s_sequences_fitting_annotation.txt' % (outdir, enzyme), 'w') as outfile:
        outfile.write('seq\n')
        for sequence in target_sequence_overlapping_exon:
            if enzyme == 'Cas9':
                if 'NGG' in sequence:
                    gRNA = sequence.replace('NGG', '')

                    if unique_dictionary.has_key(gRNA):
                        unique_dictionary[gRNA] = unique_dictionary[gRNA] + 1
                    else:
                        unique_dictionary[gRNA] = 1
                else:
                    continue
            else:
                gRNA = sequence.replace('TTTN', '')

                if unique_dictionary.has_key(gRNA):
                    unique_dictionary[gRNA] = unique_dictionary[gRNA] + 1
                else:
                    unique_dictionary[gRNA] = 1

        for key in unique_dictionary.keys():
            if unique_dictionary[key] == 1:
                outfile.write('%s\n' % key)
            else:
                continue

    sys.stdout.write('gRNA fitting annotation requirement for %s enzyme complete\n' % enzyme)

def gRNA_annotation_through_interval_tree(gRNA_file, interval_tree, cpf1):
    """annotate gRNAs according to genomic feature

    :param gRNA_file: absolute filepath to gRNA_file (field 1 gRNA target sequence, field 2 is coordinate)
    :param interval_tree: tree object, first output of create_inttree_from_file()
    :return: list and dictionary object, list of gRNA target sequences fulfilling annotation requirement: dictionary
             object with gRNA target sequence as key and coordinates as values

    """
    sequence_dictionary = {}
    target_dictionary,target_dictionary_annotation = defaultdict(list),defaultdict(list)
    target_sequence_overlapping_exon = []
    with open(gRNA_file, 'r') as infile:
        for line in infile:
            clean_line = line.lstrip().rstrip()
            parts = clean_line.split()
            target_sequence, target_coordinate = parts[0], parts[1]

            if sequence_dictionary.has_key(target_sequence):
                sequence_dictionary[target_sequence] = sequence_dictionary[target_sequence] + 1
            else:
                sequence_dictionary[target_sequence] = 1

            target_coordinate_parts = target_coordinate.split(':')
            chromosome, coordinate, strand = target_coordinate_parts[0], int(target_coordinate_parts[1]), \
                                             target_coordinate_parts[2]

            if cpf1 == 1:
                if strand == '+':
                    end_coordinate = coordinate + 19  # only for Cpf1
                    start_coordinate = end_coordinate
                elif strand == '-':
                    start_coordinate = coordinate - 19  # only for Cpf1
                    end_coordinate = start_coordinate
                else:
                    sys.stderr.write('%s is not a valid strand character\n' % strand)
                    continue

            else:
                if strand == '+':
                    end_coordinate = coordinate + 17 # only for Cas9
                    start_coordinate = end_coordinate
                elif strand == '-':
                    start_coordinate = coordinate - 17 # only for Cas9
                    end_coordinate = start_coordinate
                else:
                    sys.stderr.write('%s is not a valid strand character\n' % strand)
                    continue

            try:
                annotation_tree_grab = interval_tree.get(chromosome)
                annotation_tree_traversal = annotation_tree_grab.find(start_coordinate, end_coordinate)
                if annotation_tree_traversal:
                    target_dictionary[target_sequence].append(
                        '%s:%s-%s' % (chromosome, start_coordinate, end_coordinate))
                    target_dictionary_annotation[target_sequence].append('%s:%s-%s_%s' % (chromosome, start_coordinate, end_coordinate,annotation_tree_traversal))
                else:
                    continue

            except AttributeError:
                sys.stderr.write('%s has attribute not found in interval tree\n' % (target_coordinate))
                continue

    for key in target_dictionary.keys():
        if sequence_dictionary[key] == 1:
            target_sequence_overlapping_exon.append(key)
        else:
            continue

    sys.stdout.write('coordinate annotation query complete\n')
    return target_sequence_overlapping_exon,target_dictionary,target_dictionary_annotation

def create_inttree_from_file(infile):
    """Create interval tree to store annotations

    Args:
    infile: handle of open BED file with annotations

    Return:
    dictionary {chromosome name : interval tree with coordinates}
    """
    genome = {}
    for line in infile:
        clean_line = line.strip()
        parts = clean_line.split()
        chrom, start, stop = parts[0], int(parts[1]), int(parts[2])
        name = parts[3]
        tree = None
        #if chromosome already in tree, index to this tree
        if chrom in genome:
            tree = genome[chrom]
        else:
            #first time we encounter chromosome, create a new interval tree
            tree = IntervalTree()
            genome[chrom] = tree
        #add interval to tree
        tree.add(start, stop, name)
    return genome

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():
    """
    cpf1=1
    gRNA_file = '/Users/pereza1/Projects/Ventura/CRISPR/data/dm6/dm6_all_kmers.txt'
    annotations = '/Users/pereza1/Projects/Ventura/CRISPR/code/crispr-project/guidescan-crispr/guidescan/annotation_bed/dm6/dm6_exons_completeAnnotation.bed'
    outdir = '/Users/pereza1/Desktop'

    interval_tree = create_interval_tree(annotations)
    x = interval_tree.get('chr4')
    x.find(1050297,1053703)

    cpf1 = 1
    gRNA_file = '/Users/pereza1/Projects/Ventura/CRISPR/data/dm6/dm6_all_kmers.txt'
    annotations = '/Users/pereza1/Projects/Ventura/CRISPR/code/crispr-project/guidescan-crispr/guidescan/annotation_bed/dm6/dm6_exons_completeAnnotation.bed'
    outdir = '/Users/pereza1/Desktop'

    interval_tree = create_interval_tree(annotations)
    x = interval_tree.get('chr4')
    x.find(1047640, 1047640)
    """

    # user inputs
    annotations,gRNA_file,outdir,cpf1 = arg_parser()

    # generate interval tree
    interval_tree = create_interval_tree(annotations)

    # annotation of target sites with interval tree
    target_sequence_overlapping_exon,target_dictionary,target_dictionary_annotation = gRNA_annotation_through_interval_tree(gRNA_file,interval_tree,cpf1)

    # pickle dictionary of sequences and coordinates
    with open('%s/gRNA_target_coordinates_dict.pkl' % outdir, 'w') as out_pickle:
        pickle.dump(target_dictionary,out_pickle)

    with open('%s/gRNA_target_coordinates_annotation_dict.pkl' % outdir,'w') as out_annotation_pickle:
        pickle.dump(target_dictionary_annotation,out_annotation_pickle)

    # write out
    gRNA_w_annotations_writeout(target_sequence_overlapping_exon,outdir,cpf1)

    # end message for user
    sys.stdout.write('annotation selection of gRNAs for %s complete\n' % gRNA_file)

if __name__ == '__main__':
    main()

