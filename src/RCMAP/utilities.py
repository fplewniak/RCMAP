from RCMAP.Alignment import Alignments

def get_positions_list(positions):
    """
    :param positions:  #like ['3:10', '8:25' ,'32', '45:']
    :return: #like [[3,10], [8,25], [32], [45,None]]
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
            positions_list.append([int(positions[k])])
    return positions_list
