from Bio import AlignIO

from RCMAP.Classification_AA import AAcategories

file = "ArsM_aln.faa"  #Definition du fichier

def read_alignments(file,seqs_to_evaluate=["Q968Z2","WP_045226361"]):
    alignment = AlignIO.read(file, "fasta")
    seqeval = alignment[record.id == "Q968Z2"or "WP_045226361"]
    return seqeval

#def get_cat_at_pos(seqref,pos):
    #lire set(aa) à la position dans chaque séquence de seqref
    #return (find_category(set(aa)))

#def get_cat_in_range(seqref,debut=None,fin=None):
    #return [get_cat_at_pos(seqref,pos) for pos in range(debut,fin)]


print(read_alignments(file))