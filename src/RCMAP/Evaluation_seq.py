from RCMAP.Alignment import Alignments

file = "ArsM_aln.faa"
sequences_to_evaluate = ["WP_045226361.1", "Q969Z2"]
positions = ['3:10', '8:25' ,'32', '45:']

(seqref, seqeval) = Alignments(file,sequences_to_evaluate).get_alignments()
list_of_positions = Alignments(file,sequences_to_evaluate).get_positions_list(positions)
list_of_categories =  Alignments(file,sequences_to_evaluate).get_category_list(list_of_positions)

