from RCMAP.utilities import sort_categories


class AAcategories:
    """
    Defines the categories and determines in which category a set belongs to.
    """

    def __init__(self):
        self.categories = {'Negative': set("DE"), 'Polar': set("DEKRHQNSCTYW"), 'Tiny': set("ASCG"),
            'Positive': set("RKH"), 'Aliphatic': set("IVL"), 'Aromatic': set("FYWH"),
            'Charged': set("DEKRH"), 'Small': set("ASCGPNDTV"),'Hydrophobic': set("IVLFYWHMKTGAC")}

    def find_category(self, ens):
        """
        Finds the smallest category in which a set of amino acids (it could
        also be a unique amino acid) is included.
        Removes the gaps in the set and treats the amino acids 'B' and 'Z' which are ambiguous:
        'B' represents the amino acids 'D' or 'N', 'Z' represents the amino acids 'Q' or 'E'.

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
        for category in list_cat:
            if ens.issubset(self.categories[category]):
                return self.categories[category], category
        return set("IVLFYWHMKTGACPSNDEQR"), 'Any'

