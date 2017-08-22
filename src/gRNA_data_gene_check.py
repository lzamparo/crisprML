__author__ = 'Alexendar Perez'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""check which genes the gRNAs used for training data in CRISPR ML task are hitting in mm10"""

#################
#               #
#   Libraries   #
#               #
#################

import sys
import pickle
import argparse

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pickle',help='absolute filepath to gRNA_target_coordinates_annotation_dict.pkl',required=True)
    parser.add_argument('-i','--infile',help='absolute filepath to ultra.txt',required=True)
    parser.add_argument('-o','--outdir',help='absolute filepath to output directory',required=True)
    parser.add_argument('--cpf1',help='indication if cpf1 is being queried',default=None)

    args = parser.parse_args()
    ultra = args.infile
    in_pickle = args.pickle
    outdir = args.outdir
    cpf1 = args.cpf1

    return ultra,in_pickle,outdir,cpf1

def gene_symbol_writeout(outdir, gene_set):
    """write out file with gene symbols that are targeted

    :param outdir: absolute filepath to output directory
    :param gene_set: first output of determine_gene_targets()
    :return: writes out a file with gene symbols that are targeted

    """
    with open('%s/gene_set.txt' % (outdir), 'w') as outfile:
        for gene in gene_set:
            outfile.write('%s\n' % gene)

    sys.stdout.write('gene symbols for targeted genes written out\n')

def determine_gene_targets(ultra,gRNA_target_coordinates_annotation_dict,cpf1=None):
    """determine how many genes are targeted by ClassTask gRNAs

    :param ultra: list of gRNA sequences to be used for ClassTask
    :param gRNA_target_coordinates_annotation_dict: first output of load_pickle()
    :return: set object with gene symbols

    """
    gene_lst = []
    with open(ultra, 'r') as infile:
        for line in infile:
            clean_line = line.lstrip().rstrip()
            parts = clean_line.split()
            if cpf1:
                gRNA = 'TTTN%s' % parts[0]
            else:
                gRNA = '%sNGG' % parts[0]
            gene = gRNA_target_coordinates_annotation_dict[gRNA]
            if gene:  # coding exon
                for i in gene:
                    symbol = i.split("_['")[1].split('_')[0]
                    gene_lst.append(symbol)
            else:  # UTR
                continue

    gene_set = set(gene_lst)
    sys.stdout.write('%s genes are targeted' % (len(gene_set)))

    return gene_set

def load_pickle(in_pickle):
    """deserialize pickle

    :param in_pickle: absolute filepath to serialized pickle object
    :return: deserialized pickle object

    """
    with open(in_pickle, 'r') as p:
        data = pickle.load(p)

    sys.stdout.write('pickle object deserialized\n')
    return data

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():
    """
    outdir = '/Users/pereza1/Projects/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cas9/split/human_trie_processed/aggregate/genes_hit'
    ultra = '/Users/pereza1/Projects/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cas9/split/human_trie_processed/aggregate/filtered/ultra.txt'
    in_pickle = '/Users/pereza1/Projects/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cas9/augmented_annotations/gRNA_target_coordinates_annotation_dict.pkl'
    """

    # user inputs
    ultra, in_pickle, outdir, cpf1 = arg_parser()

    # deserialize pickle
    gRNA_target_coordinates_annotation_dict = load_pickle(in_pickle)

    # get gRNAs that meet processing requirements
    gene_set = determine_gene_targets(ultra,gRNA_target_coordinates_annotation_dict,cpf1)

    # write out value of gene set to file
    gene_symbol_writeout(outdir,gene_set)

    # user end message
    sys.stdout.write('gene intersection processing complete\n')

if __name__ == '__main__':
    main()


