class AAcategories:

    # Initializer / Instance Attributes
    def __init__(self):
        self.categories = {
            'Negative' : set("DE") ,'Polar' : set("DEKRHQNSCTYW"),'Tiny' : set("ASCG") ,
            'Positive' : set("RK") ,
            'Aliphatic' : set("IVL") ,
            'Aromatic' : set("FYWH"),
            'Charged' : set("DEKRH"),
            'Small' : set("ASCGPNDTV"),
            'Hydrophobic' : set("IVLFYWHMKTGAC")
        }
    def sort_categories(self):
        return ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar', 'Hydrophobic']

    def find_category(self,ens):
        """This function returns the smallest category in which the set is included"""
        if "-" in ens :
            ens.remove("-")
        if len(ens) <= 1 :
            return(ens)
        List = self.sort_categories()
        for k in range(len(self.categories)):
            if ens.issubset(self.categories[List[k]]):
                return(self.categories[List[k]])
        return(set("IVLFYWHMKTGACPSNDEQR"))

    def compatibility(self,AA,category):
        if AA < category:
            return True
        return False


ens = {"-","-","K","D"}
AA = {"F"}
category = set("DEKRH")
objet = AAcategories()
#print(objet.find_category(ens))
print(objet.compatibility(AA,category))