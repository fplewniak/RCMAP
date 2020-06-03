from RCMAP.Alignment import Alignments
from RCMAP.Utilities import get_positions_list, find_key
from RCMAP.Classification_AA import AAcategories

file = "ArsM_aln.faa"
sequences_to_evaluate = ["WP_045226361.1", "Q969Z2"]
positions = ['2:4','40', '50:52']

alignments = Alignments(file, sequences_to_evaluate)
list_of_positions = get_positions_list(positions)
list_of_categories = alignments.get_category_list(list_of_positions)

# Nouvelle version
for s in sequences_to_evaluate:
    print(s)
    for start, end in list_of_positions:
        for i in range(start, end):
            aa = alignments.get_aa_at_pos(i, s)
            print(i, " : ", aa, " ; Compatibility : ",
                  AAcategories().compatibility(set(aa), list_of_categories[i]), " ; ",
                  find_key(list_of_categories[i]), alignments.list_of_aa_ref[i-1])

