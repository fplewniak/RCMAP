def get_positions_list(positions):
    """
    :param positions:  positions like ['3:10', '8:25' ,'32', '45:']
    :return: positions like [[3,10], [8,25], [32], [45,None]]
    """
    positions_list = []
    for k in range(len(positions)):
        l = []
        if ':' in positions[k]:
            r = positions[k].rpartition(':')
            if r[0] != '':
                l.append(int(r[0]))
            else:
                l.append(None)
            if r[2] != '':
                l.append(int(r[2]))
            else:
                l.append(None)
            positions_list.append(l)
        else:
            positions_list.append([int(positions[k]), int(positions[k])])
    return positions_list


def summary_info(List_compatibility, List_info):
    count_False, count_True = 0, 0
    info_False, info_True = 0, 0
    for k in range(len(List_compatibility)):
        if List_compatibility[k] is True:
            count_True += 1
            info_True += List_info[k]
        else:
            count_False += 1
            info_False += List_info[k]
    return count_True, round(info_True, 2), count_False, round(info_False, 2)
