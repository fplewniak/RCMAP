
def sort_categories():
    """
    :return: the list of the categories
    """
    return ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar',
            'Hydrophobic']


def get_positions_list(positions):
    """
    :param positions:  positions like ['3:10', '8:25' ,'32', '45:']
    :return: positions like [[3,10], [8,25], [32], [45,None]]
    """
    positions_list = []
    for pos in positions:
        new_pos = []
        if ':' in pos:
            rep = pos.rpartition(':')
            if rep[0] != '':
                new_pos.append(int(rep[0]))
            else:
                new_pos.append(None)
            if rep[2] != '':
                new_pos.append(int(rep[2]))
            else:
                new_pos.append(None)
            positions_list.append(new_pos)
        else:
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
