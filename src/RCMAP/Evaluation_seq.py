from RCMAP.Alignment import Alignments
from RCMAP.Utilities import get_positions_list, find_key
from RCMAP.Classification_AA import AAcategories

file = "ArsM_aln.faa"
sequences_to_evaluate = ["WP_045226361.1", "Q969Z2"]
positions = ['2:4','40', '50:52']

alignments = Alignments(file, sequences_to_evaluate)
list_of_positions = get_positions_list(positions)
list_of_categories = alignments.get_category_list(list_of_positions)

for s in range(len(sequences_to_evaluate)):
    list_of_aa_in_seqeval = alignments.get_aa_list(list_of_positions, sequences_to_evaluate[s])
    print(sequences_to_evaluate[s])
    for k in range(len(list_of_categories)):
        pos = list_of_positions[k][0]
        for i in range(len(list_of_categories[k])):
            print(pos, " : ", list_of_aa_in_seqeval[k][i], " ; Compatibility : ",
                  AAcategories().compatibility(list_of_aa_in_seqeval[k][i],
                                               list_of_categories[k][i])
                  , " ; ",find_key(list_of_categories[k][i]), alignments.list_of_aa_ref[pos-1])
            pos += 1

# Nouvelle version
for s in sequences_to_evaluate:
    print(s)
    for k in range(len(list_of_positions)):
        category = alignments.list_of_categories[
                   list_of_positions[k][0] - 1:list_of_positions[k][1]]
        pos = 0
        for i in range(list_of_positions[k][0], list_of_positions[k][1] + 1):
            aa = alignments.get_aa_at_pos(i, s)
            print(i, " : ", aa, " ; Compatibility : ",
                  AAcategories().compatibility(set(aa), category[pos]), " ; ",
                  find_key(category[pos]), alignments.list_of_aa_ref[i-1])
            pos += 1
