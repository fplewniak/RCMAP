from Bio import AlignIO

alignment = AlignIO.read("data/ArsM_aln.faa", "fasta")
print(alignment)

for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))

