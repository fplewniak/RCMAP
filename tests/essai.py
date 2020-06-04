def get_preserved_area(list_pos):
    """
    :param list_pos: list of the positions where the information is sufficiente
    :return: the preserved areas
    """
    list_conserved_pos = []
    k = 0
    while k < len(list_pos) - 1:
        if list_pos[k] == list_pos[k + 1] - 1:
            i = k
            while k < len(list_pos) - 1 and list_pos[k] == list_pos[k + 1] - 1:
                k += 1
            list_conserved_pos.append([list_pos[i], list_pos[k]])
        k += 1
    return list_conserved_pos


list_pos = [2, 4, 5, 20, 28, 30, 31, 32, 33, 80]
print(get_preserved_area(list_pos))


def frequence_aa(self, aa, pos):
    if self.count_aa_ref()[pos - 1][aa] == 0:
        return 0
    return self.count_aa_ref()[pos - 1][aa] / len(self.seqrefs)