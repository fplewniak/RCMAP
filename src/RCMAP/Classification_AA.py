class AAcategories:

    # Initializer / Instance Attributes
    def __init__(self):
        self.categories = {
            'Negative' : set(["D","E"]) ,
            'Positive' : set(["R","K"]) ,
            'Aliphatic' : set(["I","V","L"]) ,
            'Tiny' : set(["A", "S","C","G"]) ,
            'Aromatic' : set(["F", "Y", "W", "H"]),
            'Charged' : set(["D", "E", "K", "R","H"]),
            'Small' : set(["A", "S","C","G","P","N","D","T","V"]),
            'Polar' : set(["D", "E", "K", "R","H","Q","N","S","C","T","Y","W"]),
            'Hydrophobic' : set(["I","V","L","F", "Y", "W", "H","M","K","T","G","A","C"])
        }

    def find_category(self,ens):
        """This function returns the smallest category in which the set is included"""
        if "-" in ens :
            ens.remove("-")
        if len(ens) == 1 :
            return(print(ens))
        List = ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar','Hydrophobic']
        for k in range(len(self.categories)):
            if ens.issubset(self.categories[List[k]]):
                return(print(List[k]))
        return(print(set(["I","V","L","F", "Y", "W", "H","M","K","T","G","A","C","P","S","N","D","E","Q","R"])))


ens = {"G","-","C"}
objet = AAcategories()
objet.find_category(ens)

