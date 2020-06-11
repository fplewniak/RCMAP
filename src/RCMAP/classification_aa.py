from RCMAP.utilities import sort_categories


class AAcategories:
    """
    AAcategories is used to define the categories and test if an amino acid belongs to a category
    """

    # Initializer / Instance Attributes
    def __init__(self):
        self.categories = {
            'Negative': set("DE"), 'Polar': set("DEKRHQNSCTYW"), 'Tiny': set("ASCG"),
            'Positive': set("RKH"),
            'Aliphatic': set("IVL"),
            'Aromatic': set("FYWH"),
            'Charged': set("DEKRH"),
            'Small': set("ASCGPNDTV"),
            'Hydrophobic': set("IVLFYWHMKTGAC")
        }

    def find_category(self, ens):
        """
        :param ens: set
        :return: the smallest category in which the set is included, and its name
        """
        if "-" in ens:
            ens.remove("-")
        if "B" in ens:
            ens.remove("B")
            ens.add("D")
            ens.add("N")
        if "Z" in ens:
            ens.remove("Z")
            ens.add("Q")
            ens.add("E")
        if len(ens) <= 1:
            return ens, ens
        list_cat = sort_categories()
        for k in range(len(self.categories)):
            if ens.issubset(self.categories[list_cat[k]]):
                return self.categories[list_cat[k]], list_cat[k]
        return set("IVLFYWHMKTGACPSNDEQR"), 'Any'

