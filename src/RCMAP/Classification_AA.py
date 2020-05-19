class AAcategories:

    # Initializer / Instance Attributes
    def __init__(self):
        self.categories = {
            'Negative' : set(["D","E"]) ,
            'Positive' : set(["R","K"]) ,
            'Aliphatic' : set(["I","V","L"]) ,
            'Tiny' : set(["A", "S","Csh","G"]) ,
            'Aromatic' : set(["F", "Y", "W", "H"]),
            'Charged' : set(["D", "E", "K", "R","H"]),
            'Small' : set(["A", "S","Csh","G","P","N","D","T","V","Css"]),
            'Polar' : set(["D", "E", "K", "R","H","Q","N","S","Csh","T","Y","W"]),
            'Hydrophobic' : set(["I","V","L","F", "Y", "W", "H","M","K","T","Cs","G","A","Csh"])
        }

    def find_category(self,ens):
        """This function returns the smallest category in which the set is included"""
        if len(ens) == 1 :
            return(print("The class contains only one amino acid : ",list(ens)[0]))
        List = ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar','Hydrophobic']
        for k in range(len(self.categories)):
            if ens.issubset(self.categories[List[k]]):
                return(print("The set is contained in the category : ",List[k]))
        return(print("The set is not contained in any of the categories"))


ens = {"P"}
objet = AAcategories()
objet.find_category(ens)


