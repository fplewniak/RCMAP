import numpy

def sort_categories():
    """
    :return: the list of the categories
    """
    return ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar',
            'Hydrophobic']


def compatibility(amino_acid, category, gaps=False):
    """
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


def get_positions_list(positions, posmax):
    """
    :param posmax: maximal authorized value for a position (length of the sequence)
    :param positions:  positions like ['3:10', '8:25' ,'32', '45:', ':5']
    :return: positions like [[3,10], [8,25], [32,32], [45,posmax], [1:5]]
    """
    positions_list = []
    for pos in positions:
        new_pos = []
        if ':' in pos:
            rep = pos.rpartition(':')
            if rep[0] != '':
                new_pos.append(int(rep[0]))
            else:
                new_pos.append(1)
            if rep[2] != '':
                new_pos.append(int(rep[2]))
            else:
                new_pos.append(posmax)
            if new_pos[0] > new_pos[1]:
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
    * The Bartlett window is defined as : w(n) = 2/(M-1) * ((M-1)/2 - abs(n-(M-1)/2))
    * The Hamming window is defined as : w(n) = 0.54 - 0.46 * cos((2*pi*n)/(M-1))  0 <= n <= M-1
    * The Hanning window is defined as : w(n) = 0.5 - 0.5 * cos((2*pi*n)/(M-1))  0 <= n <= M-1
    * The flat window is defined as : w(n) = 1
    :param window_method: Calculation method of the weights at every position in the window
    :return: the list of weights
    """
    if window_method == 'flat':
        w = numpy.ones(window)
    else:
        w = eval('numpy.' + window_method + '(window)')
    return w


