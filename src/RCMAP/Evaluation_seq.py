from RCMAP.Alignment import Alignments
from RCMAP.Utilities import get_positions_list
from RCMAP.Classification_AA import AAcategories

file = "ArsM_aln.faa"
sequences_to_evaluate = ["WP_045226361.1", "Q969Z2"]
positions = ['2:4','40','50:52']

(seqref, seqeval) = Alignments(file,sequences_to_evaluate).get_alignments()
list_of_positions = get_positions_list(positions)
list_of_categories = Alignments(file,sequences_to_evaluate).get_category_list(list_of_positions)
print(list_of_positions)
print(list_of_categories)
list_compatibility_all = []

for s in seqeval:
    name_seqeval = s.id
    list_of_aa_in_seqeval = Alignments(file,sequences_to_evaluate).get_aa_list_in_seqeval(name_seqeval,list_of_positions)
    print(name_seqeval,list_of_aa_in_seqeval)
    list_compatibility = []
    for k in range(len(list_of_categories)):
        pos = list_of_positions[k][0]
        for i in range(len(list_of_categories[k])):
            list_compatibility.append(pos)
            list_compatibility.append(AAcategories().compatibility(list_of_aa_in_seqeval[k][i], list_of_categories[k][i]))
            pos+=1
    list_compatibility_all.append([name_seqeval,list_compatibility])

print(list_compatibility_all)
