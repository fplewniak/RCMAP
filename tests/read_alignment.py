from Bio import AlignIO


class Read_sequences:

    def __init__(self,n=6):
        """Save the number of reference sequences"""
        self.n= n

    def reference_sequences(self):
        self.alignment = AlignIO.read("ArsM_aln.faa", "fasta")
        self.ref_sequences = self.alignment[0:self.n]
        return (self.ref_sequences)
        return len(self.ref_sequences[0])

    def create_set(self):
        self.List_set =[]
        for i in range (len(self.ref_sequences[0])):
            set_AA = set()
            for k in range (self.n):
                set_AA.add(self.alignment[k][i])
            self.List_set.append(set_AA)
        return self.List_set





objet = Read_sequences()
print(objet.reference_sequences())
print(objet.create_set())