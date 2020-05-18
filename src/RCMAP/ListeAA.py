Categories = {'Proline' : set(["P"]) ,
            'Negative' : set(["D","E"]) ,
            'Positive' : set(["R","K"]) ,
            'Aliphatic' : set(["I","V","L"]) ,
            'Tiny' : set(["A", "S","Csh","G"]) ,
            'Aromatic' : set(["F", "Y", "W", "H"]),
            'Charged' : set(["D", "E", "K", "R","H"]),
            'Small' : set(["A", "S","Csh","G","P","N","D","T","V","Css"]),
            'Polar' : set(["D", "E", "K", "R","H","Q","N","S","Csh","T","Y","W"]),
            'Hydrophobic' : set(["I","V","L","F", "Y", "W", "H","M","K","T","Cs","G","A","Csh"])}

#print(Categories)
#print(Categories['Negative'])
#print(Categories['Charged'])
#print(Categories['Negative'].issubset(Categories['Charged']))

ens = {"V","I","Q"}
#print(len(Categories))

def find_category(ens):
    """Cette fonction renvoie tous les ensembles dans lesquels ens est inclu"""
    n=len(Categories)
    list = ['Proline','Negative','Positive','Aliphatic','Tiny','Aromatic','Charged','Small','Polar','Hydrophobic']
    for k in range (n) :
        if ens.issubset(Categories[list[k]]) :
            print (list[k])

def find_category1(ens):
    """Cette fonction renvoie le plus petit ensemble dans lequel ens est inclu"""
    n=len(Categories)
    i=0
    list = ['Proline','Negative','Positive','Aliphatic','Tiny','Aromatic','Charged','Small','Polar','Hydrophobic']
    for k in range (n) :
        if ens.issubset(Categories[list[k]]) :
            print (list[k])
            i=1
            break
    if i==0:
        print("L'ensemble n'est contenu dans aucune des catégories")

#find_category(ens)
find_category1(ens)