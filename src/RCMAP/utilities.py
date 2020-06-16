import numpy


def sort_categories():
    """
    Function used in classification_aa.py in find_category() which needs the categories sorted
    by size.
    :return: the list of sorted categories, from that which contains the least of amino acids,
    to that which contains the most
    """
    return ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar',
            'Hydrophobic']


def compatibility(amino_acid, category, gaps=False):
    """
    This function is used to test if a set of amino acid (or a single amino acid) is included in a
    category. It also treat the amino acids 'B' and 'Z' which are ambiguous. They represent two
    possibilities of amino acids.
    If there is a gap in the set to test, the function returns gaps (True or False). The parameter
    gaps is defined in the main function in evaluation_seq.py. When compatibility() is used,
    gaps = True when the parameter gaps is defined True by the user AND there is a gap in the set
    of amino acids observed in the reference sequences.
    :param gaps: consider gaps if gaps is True
    :param amino_acid: set of amino acids
    :param category: a category (set)
    :return: True if the set of amino acids is included in the category, False if not
    """
    if "B" in amino_acid:
        amino_acid = {"D", "N"}
    if "Z" in amino_acid:
        amino_acid = {"Q", "E"}
    if "-" in amino_acid:
        return gaps
    if amino_acid <= category:
        return True
    return False


def get_positions_list(positions, pos_max):
    """
    This function is used to transform the positions given by the user in a usable form for the
    other functions. It replaces non given positions by the beginning or the end of the sequence.
    :param pos_max: maximal authorized value for a position (length of the sequence)
    :param positions:  positions like ['3:10', '8:25' ,'32', '45:', ':5', ':']
    :return: positions like [[3,10], [8,25], [32,32], [45,pos_max], [1,5], [1,pos_max]]
    """
    positions_list = []
    for pos in positions:
        new_pos = []
        if ':' in pos:  # The position is a interval, not a single position
            rep = pos.rpartition(':')
            if rep[0] != '':    # If the beginning of the interval is given
                new_pos.append(int(rep[0]))
            else:
                new_pos.append(1)
            if rep[2] != '':    # If the end of the interval is given
                new_pos.append(int(rep[2]))
            else:
                new_pos.append(pos_max)
            if new_pos[0] > new_pos[1]:   # Reset positions from smallest to largest
                new_pos[0], new_pos[1] = new_pos[1], new_pos[0]
            if new_pos[0] <= 0:
                print('Error in positions, negative value : {pos}'.format(pos=new_pos[0]))
                exit(1)
            positions_list.append(new_pos)
        else:
            if int(pos) <= 0:
                print('Error in positions, negative value : {pos}'.format(pos=int(pos)))
                exit(1)
            positions_list.append([int(pos), int(pos)])
    return positions_list


def summary_info(list_compatibility, list_info):
    """
    This function is used to resume the global information of a sequence. It resume the total
    number of 'True' positions (which means the amino acid of the evaluated sequence is compatible
    with the category observed in the reference sequences) and the total of information associated.
    In the same way for 'False' positions.
    :param list_compatibility: list of compatibility True or False at every position
    :param list_info: list of information at every position
    :return: number of True, sum of information for True positions, number of False, sum of
    information for False positions
    """
    count_false, count_true = 0, 0
    info_false, info_true = 0, 0
    for k in range(len(list_compatibility)):
        if list_compatibility[k] is True:
            count_true += 1
            info_true += list_info[k]
        else:
            count_false += 1
            info_false += list_info[k]
    return count_true, round(info_true, 2), count_false, round(info_false, 2)


def get_weight(window, window_method):
    """
    This function aims to define the weights of positions in a window on the information of a
    position whose information is evaluated. These weights can be calculated with different method:
    * The Bartlett window is defined as : w(n) = 2/(M-1) * ((M-1)/2 - abs(n-(M-1)/2))
    * The Hamming window is defined as : w(n) = 0.54 - 0.46 * cos((2*pi*n)/(M-1))  0 <= n <= M-1
    * The Hanning window is defined as : w(n) = 0.5 - 0.5 * cos((2*pi*n)/(M-1))  0 <= n <= M-1
    * The flat window is defined as : w(n) = 1
    :param window: length of the window
    :param window_method: Calculation method of the weights at every position in the window
    :return: the list of weights
    """
    if window_method == 'flat':
        w = numpy.ones(window)
    else:
        w = eval('numpy.' + window_method + '(window)')
    return w
