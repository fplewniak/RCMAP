import argparse
from RCMAP.alignment import Alignments
from RCMAP.utilities import get_positions_list, summary_info, compatibility


def get_params():
    """
    Get parameters from command line :
    * file : the file which contains the alignments of the reference sequences and the sequences to
      evaluate
    * seqeval : the name of the sequences to evaluate
    * positions : the positions or interval on which the analysis must be made, begins at 1
    * method : the calculation method of the frequencies of the amino acids in the background
      entropy
    * gaps : True or false. It depends if the user want to consider gaps. If gaps = True, when there
      is a gap in the evaluated sequence and a gap in the amino acids observed in the reference
      sequences, the compatibility is evaluated as True. Otherwise not
    * strict : True or False. True if you want to compare the amino acid in the evaluated sequence
      only with the amino acids observed in the reference sequences and not with the associated
      category
    * min_info : Minimum of information required by the user for display
    * window : Number of positions around the evaluated position to calculate the average of its
      information. It must be odd. Used to consider the environment of a position to calculate the
      information carried
    * window_method : Calculation method of the weights of the positions around the position which
    information is evaluated using a window

    return: parameters
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--file', metavar='File', help='Alignment file', required=True)
    parser.add_argument('--seqeval', metavar='String', nargs='+', help='Sequences to evaluate',
                        required=True)
    parser.add_argument('--positions', nargs='+', help='List of positions to evaluate',
                        default=[':'])
    parser.add_argument('--method', type=str,
                        choices=['equiprobable', 'database', 'ref'], help='Calculation method of '
                                                                          'the background entropy'
                                                                          ' for the information',
                        default='equiprobable')
    parser.add_argument('--gaps', help='True if you want to consider the gaps', action='store_true')
    parser.add_argument('--strict', help='True if you want to compare the amino acid in the '
                                         'evaluated sequence only with the amino acids observed '
                                         'in the reference sequences', action='store_true')
    parser.add_argument('--min_info', type=float,
                        help='Minimum of information required for display', default=0)
    parser.add_argument('--window', type=int,
                        help='Number of positions to calculate the average of '
                             'information, must be odd', default=1)
    parser.add_argument('--window_method', type=str,
                        choices=['bartlett', 'hamming', 'hanning', 'flat'],
                        help='Calculation method of the weights of positions to calculate the '
                             'information of a position, using a window',
                        default='flat')

    params = parser.parse_args()
    if params.window % 2 == 0:
        print('Window must be odd, {window} is not'.format(window=params.window))
        exit(1)
    if params.window <= 0:
        print('Window must be positive, {window} is not'.format(window=params.window))
        exit(1)
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
    list_of_positions = get_positions_list(params.positions, len(alignments.seqeval[0]))
    for seq in params.seqeval:
        print('\n')
        list_compatibility, list_info = [], []
        for start, end in list_of_positions:
            for i in range(start, end + 1):
                info = round(
                    alignments.information_pos(i, params.method, params.gaps, params.window,
                                               params.window_method), 2)
                if info >= params.min_info:
                    print('{name_seq} : {pos} : {aa} : {test} : {info:.2f} : {cat} {obs}'.format(
                        name_seq=seq.rjust(0), pos=str(i).rjust(4),
                        aa=alignments.get_aa_at_pos(i, seq),
                        test=str(compatibility(set(alignments.get_aa_at_pos(i, seq)),
                                               alignments.get_cat_at_pos(i, params.strict),
                                               params.gaps & ('-' in {aa for aa in
                                                                      alignments.aa_observed[
                                                                          i - 1]}))).rjust(6),
                        cat=str(alignments.get_cat_name_at_pos(i)).center(10),
                        obs=alignments.get_aa_observed_at_pos(i),
                        info=info))

                    list_compatibility.append(
                        compatibility(set(alignments.get_aa_at_pos(i, seq)),
                                      alignments.get_cat_at_pos(i, params.strict),
                                      params.gaps & ('-' in {aa for aa in
                                                             alignments.aa_observed[i - 1]})))
                    list_info.append(
                        alignments.information_pos(i, params.method, params.gaps, params.window,
                                                   params.window_method))

        print('\n', seq, '\n', 'Number of True '.rjust(0),
              summary_info(list_compatibility, list_info)[0],
              ' : Information True ', summary_info(list_compatibility, list_info)[1], '\n',
              'Number of False '.rjust(
                  0), summary_info(list_compatibility, list_info)[2],
              ' : Information False ', summary_info(list_compatibility, list_info)[3])


if __name__ == '__main__':
    main()
