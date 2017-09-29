__author__ = 'Alexendar Perez'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""extract candidate gRNAs for cutting efficiency screen"""

#################
#               #
#   Libraries   #
#               #
#################

import sys
import argparse
import pdb

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile',help='absolute filepath to output of specificity_score_distance_neighbors.py file',required=True)
    parser.add_argument('-o','--outdir',help='absolute filepath to output directory',required=True)
    parser.add_argument('--enzyme',help='enter Cas9 or Cpf1, default = Cas9',default='Cas9')

    args = parser.parse_args()
    in_file = args.infile
    outdir = args.outdir
    enzyme = args.enzyme

    return in_file,outdir,enzyme

def string_set(string_file):
    """get set of strings

    :param string_file: absolute filepath to single column string file
    :return: set object of strings

    """
    string_lst = []
    with open(string_file, 'r') as infile:
        for line in infile:
            clean_line = line.lstrip().rstrip()
            parts = clean_line.split()
            string_lst.append(parts[0])

    string_set_return = set(string_lst)
    return string_set_return

def get_gRNAs(db):
    """get gRNAs

    :param db: ultra or g_ultra list
    :return: set object with sequences as elements

    """
    hold = []
    for j in range(len(db)):
        hold.append(db[j].split(',')[0])
    hold_set = set(hold)

    return hold_set

def extract_candidate_gRNAs(in_file,enzyme):
    """extract candidate gRNAs for cutting efficiency screen

    :param in_file: absolute filepath to input file
    :return: list object, gRNAs with no near matches <= 3 starts with G,
                          gRNAs with no near matches <= 3,
                          gRNAs with no near matches <= 2
                          gRNAs with no near matches <= 1
                          gRNAs with no perfect matches

    """
    lst_c, lst_b, lst_a, ultra, g_ultra = [], [], [], [], []
    with open(in_file, 'r') as infile:
        for line in infile:
            clean_line = line.lstrip().rstrip()
            parts = clean_line.split(',')

            if enzyme == 'Cas9':
                enz_val = 'inf'
            elif enzyme == 'Cpf1':
                enz_val = '0.0'
            else:
                sys.stderr.write('%s not recognized: enter either Cas9 or Cpf1' % enzyme)
                return

            try:

                if parts[1].strip() == enz_val: #Cas9 = inf, cpf1 = 0.0
                    if float(parts[2]) == 0.0 and float(parts[3]) == 0.0 and float(parts[4]) == 0.0 and float(
                            parts[5]) == 0.0:
                        if parts[0][0] == 'G':
                            g_ultra.append(line)
                            ultra.append(line)
                        else:
                            ultra.append(line)  # no duplicate or near neighbors within Hamming distance 3
                    else:
                        pass
                else:
                    pass

            except ValueError:
                sys.stderr.write('skipping %s\n' % line)
                continue

    sys.stdout.write('gRNA extraction for %s complete\n' % (in_file))
    return g_ultra, ultra, lst_a, lst_b, lst_c

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():
    """
    in_file = '/Users/pereza1/Projects/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cas9/GuideScan_database/split/human_trie_processed/aggregate/raw_features_computed_cas9_guidescan_gRNAs_aggreagate.csv'
    outdir = '/cbio/cllab/home/aperez/Project/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cpf1/GuideScan_database/split/human_trie_processed/aggregate/filtered'
    """

    # user inputs
    in_file,outdir,enzyme = arg_parser()

    # extract candidate gRNAs
    g_ultra, ultra, lst_a, lst_b, lst_c = extract_candidate_gRNAs(in_file,enzyme)

    pdb.set_trace()
    # write out
    with open('%s/g_ultra.txt' % outdir, 'w') as g_ultra_writeout:
        for i in g_ultra:
            g_ultra_writeout.write('%s\n' % i.split()[0].strip(','))
    sys.stdout.write('%s written\n' % '%s/g_ultra.txt' % outdir)

    with open('%s/ultra.txt' % outdir, 'w') as ultra_writeout:
        for i in ultra:
            ultra_writeout.write('%s\n' % i.split()[0].strip(','))
    sys.stdout.write('%s written\n' % '%s/ultra.txt' % outdir)

    # user end message
    sys.stdout.write('gRNA extraction complete\n')

if __name__ == '__main__':
    main()

