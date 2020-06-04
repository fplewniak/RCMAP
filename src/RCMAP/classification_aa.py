class AAcategories:

    # Initializer / Instance Attributes
    def __init__(self):
        self.categories = {
            'Negative' : set("DE") ,'Polar' : set("DEKRHQNSCTYW"),'Tiny' : set("ASCG") ,
            'Positive' : set("RKH") ,
            'Aliphatic' : set("IVL") ,
            'Aromatic' : set("FYWH"),
            'Charged' : set("DEKRH"),
            'Small' : set("ASCGPNDTV"),
            'Hydrophobic' : set("IVLFYWHMKTGAC")
        }

    def sort_categories(self):
        return ['Negative', 'Positive', 'Aliphatic', 'Tiny', 'Aromatic', 'Charged', 'Small', 'Polar', 'Hydrophobic']

    def find_category(self,ens):
        """This function returns the smallest category in which the set is included, and its name"""
        if "-" in ens :
            ens.remove("-")
        if "B" in ens :
            ens.remove("B")
            ens.add("D")
            ens.add("N")
        if "Z" in ens :
            ens.remove("Z")
            ens.add("Q")
            ens.add("E")
        if len(ens) <= 1 :
            return ens, ens
        List = self.sort_categories()
        for k in range(len(self.categories)):
            if ens.issubset(self.categories[List[k]]):
                return self.categories[List[k]], List[k]
        return set("IVLFYWHMKTGACPSNDEQR"), 'Any'

    def compatibility(self,AA,category):
        if "B" in AA :
            AA = {"D","N"}
        if "Z" in AA :
            AA = {"Q","E"}
        if "-" in AA :
            return False
        if AA <= category:
            return True
        return False


