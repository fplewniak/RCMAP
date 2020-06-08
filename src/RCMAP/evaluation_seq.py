import argparse
import sys
from RCMAP.alignment import Alignments
from RCMAP.utilities import get_positions_list
from RCMAP.classification_aa import AAcategories


def get_params(argv):
    """
    Get parameters from command line
    :param argv:
    :return:
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--file', metavar='File', help='Alignment file', required=True)
    parser.add_argument('--seqeval', metavar='string', nargs='+', help='Sequences to evaluate',
                        required=True)
    parser.add_argument('--positions', nargs='+', help='List of positions to evaluate',
                        default=[':'])
    parser.add_argument('--method', type=str,
                        choices=['frq_equiprobable', 'frq_database', 'frq_alignment_ref'], default=['frq_equiprobable'])
    parser.add_argument('--gaps', help='True if you want to consider the gaps', action='store_true')
    a = parser.parse_args()
    return a


def main():
    a = get_params(sys.argv[1:])
    method = a.method
    gaps = a.gaps
    alignments = Alignments(a.file, a.seqeval)
    list_of_positions = get_positions_list(a.positions)
    list_of_categories = alignments.get_category_list(list_of_positions)
    for s in a.seqeval:
        print('\n')
        for start, end in list_of_positions:
            for i in range(1 if start is None else start,
                           len(alignments.seqeval[0]) if end is None else end + 1):
                print('{name_seq} : {pos} : {aa} : {test} : {cat} {obs} : {info}'.format(name_seq=s, pos=i,
                              aa=alignments.get_aa_at_pos(i, s),
                              test=AAcategories().compatibility(set(alignments.get_aa_at_pos(i, s)),
                                                                alignments.get_cat_at_pos(i)),
                              cat=alignments.get_cat_set_at_pos(i),
                              obs=alignments.get_aa_observed_at_pos(i), info=alignments.information_pos(i, method, gaps)))

if __name__ == '__main__':
    main()
