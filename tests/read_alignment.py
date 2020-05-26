from Bio import AlignIO

from RCMAP.Classification_AA import AAcategories

class Read_sequences:

    def __init__(self,n=6, data = "data/ArsM_aln.faa" ):
        """Save the number of reference sequences"""
        self.n= n
        self.data = data

    def reference_sequences(self):
        self.alignment = AlignIO.read(data, "fasta")
        self.ref_sequences = self.alignment[0:self.n]
        return self.ref_sequences

    def create_set_AA(self):
        self.List_AA =[]
        for i in range (len(self.ref_sequences[0])):
            set_AA = set()
            for k in range (self.n):
                set_AA.add(self.alignment[k][i])
            self.List_AA.append(set_AA)
        return self.List_AA

    def create_set_categories(self):
        self.List_set = []
        for k in range (len(self.List_AA)):
            self.List_set.append(AAcategories().find_category(self.List_AA[k]))
        return self.List_set

objet = Read_sequences()
#print(objet.reference_sequences())
#print(objet.create_set_AA())
#print(objet.create_set_categories())