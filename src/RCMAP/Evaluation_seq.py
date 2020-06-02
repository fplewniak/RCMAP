from RCMAP.Alignment import Alignments
from RCMAP.Utilities import get_positions_list
from RCMAP.Classification_AA import AAcategories

file = "ArsM_aln.faa"
sequences_to_evaluate = ["WP_045226361.1", "Q969Z2"]
positions = ['2:4', '40', '50:52']

alignments = Alignments(file, sequences_to_evaluate)
(seqref, seqeval) = alignments.get_alignments()
list_of_positions = get_positions_list(positions)
list_of_categories = alignments.get_category_list(list_of_positions)


for s in seqeval:
    list_of_aa_in_seqeval = alignments.get_aa_list(list_of_positions, s.id)
    print(s.id)
    for k in range(len(list_of_categories)):
        pos = list_of_positions[k][0]
        for i in range(len(list_of_categories[k])):
            print(pos, " : ", "in seq_eval : ", list_of_aa_in_seqeval[k][i], " ; Compatibility : ",
                  AAcategories().compatibility(list_of_aa_in_seqeval[k][i],
                                               list_of_categories[k][i])
                  , " ; category in seq_ref : ", list_of_categories[k][i])
            pos += 1

