import argparse
from RCMAP.alignment import Alignments
from RCMAP.utilities import get_positions_list, summary_info, compatibility


def get_params():
    """
    Get parameters from command line
    return: parameters
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--file', metavar='File', help='Alignment file', required=True)
    parser.add_argument('--seqeval', metavar='string', nargs='+', help='Sequences to evaluate',
                        required=True)
    parser.add_argument('--positions', nargs='+', help='List of positions to evaluate',
                        default=[':'])
    parser.add_argument('--method', type=str,
                        choices=['equiprobable', 'database', 'ref'], help='Calculation method of '
                                                                          'the background entropy'
                                                                          ' for the information',
                        default='equiprobable')
    parser.add_argument('--gaps', help='True if you want to consider the gaps', action='store_true')
    params = parser.parse_args()
    return params


def main():
    """
    :return: for each evaluated sequence :
    name of the evaluated sequence, position, amino acid in the evaluated sequence,
    compatibility of the amino acid with the category observed in reference sequences,
    information carried by the position,  the category observed in reference sequences,
    the amino acids really observed in the reference sequences.
    Summary of the information.
    """
    params = get_params()
    alignments = Alignments(params.file, params.seqeval)
    list_of_positions = get_positions_list(params.positions)
    for seq in params.seqeval:
        print('\n')
        list_compatibility, list_info = [], []
        for start, end in list_of_positions:
            for i in range(1 if start is None else start,
                           len(alignments.seqeval[0]) if end is None else end + 1):
                print('{name_seq} : {pos} : {aa} : {test} : {info:.2f} : {cat} {obs}'.format(
                    name_seq=seq.rjust(0), pos=str(i).rjust(4),
                    aa=alignments.get_aa_at_pos(i, seq),
                    test=str(compatibility(set(alignments.get_aa_at_pos(i, seq)),
                                           alignments.get_cat_at_pos(i))).rjust(6),
                    cat=str(alignments.get_cat_set_at_pos(i)).center(10),
                    obs=alignments.get_aa_observed_at_pos(i),
                    info=round(alignments.information_pos(i, params.method, params.gaps), 2)))

                list_compatibility.append(
                    compatibility(set(alignments.get_aa_at_pos(i, seq)),
                                  alignments.get_cat_at_pos(i)))
                list_info.append(alignments.information_pos(i, params.method, params.gaps))

        print('\n', seq, '\n', 'Number of True '.rjust(0),
              summary_info(list_compatibility, list_info)[0],
              ' : Information True ', summary_info(list_compatibility, list_info)[1], '\n',
              'Number of False '.rjust(
                  0), summary_info(list_compatibility, list_info)[2],
              ' : Information False ', summary_info(list_compatibility, list_info)[3])


if __name__ == '__main__':
    main()
