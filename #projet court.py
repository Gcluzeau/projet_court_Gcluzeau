#projet court

import requests
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

"""
création de la class permettant de faire fonctionner le DSSP classique
"""

class dssp_originel:
    def __init__(self, file_path):
        self.file_path = file_path

    def read(self):

        """Read PDB file and return structure."""

        parser = PDBParser(QUIET=True)

        structure = parser.get_structure('protein_structure', self.file_path)

        dssp_predictions = []

        for model in structure:

            for chain in model:

                dssp_obj = DSSP(model, self.file_path)

                for res in dssp_obj:

                    res_id, dssp_code = res[0], res[2]

                    dssp_predictions.append(dssp_code)

        return dssp_predictions


"""
Main Code
"""

"""Ouverture du ichier HB2 et selection des liaison MAIN/MAIN"""

matching_lines = []
with open("8zap.hb2", 'r') as source_file:
    for line in source_file:
        if line[33:35] == 'MM':
            matching_lines.append(line.strip())

print(len(matching_lines))

"""Création d'un dataframe avec toute les données"""

donne = []

for ligne in range(len(matching_lines)): # iteration sur element de la liste pour separer les information grace à .split()
    ligne1 = matching_lines[ligne]
    ligne1 = ligne1.split()
    donne.append(ligne1)

colonne = ['Residu1', 'Atome1', 'Residu2', 'Atome2', 'Distance', 'Type', 
    'Valeur1', 'Valeur2', 'Valeur3', 'Valeur4', 'Valeur5', 'Valeur6', 'Index']
donnes = pd.DataFrame(donne,columns = colonne)


"""Premier test d'exclusion sur l'angle """

donnes["Valeur3"] = donnes["Valeur3"].astype(float) # transforme toute les valeurs de foat réaliser des test

masque = donnes["Valeur3"] > 60 # Création du masque pour limiter l'angle de la liaison à 60°

df_filtre = donnes[masque] # Filtrer le DataFrame en utilisant le masque


"""Deuxième test d'exclusion sur la distance """


df_filtre["Distance"] = df_filtre["Distance"].astype(float) # transforme toute les valeurs de foat réaliser des test

masque2 = 0 < df_filtre["Distance"] # Création du masque pour limiter la distance et exclure le valeur négative

df_filtre2 = df_filtre[masque2]# Filtrer le DataFrame en utilisant le masque



"""Séparation des information de la chaine et du numeros de résidu"""

df_filtre2["id_res_doneur"] = df_filtre2['Residu1'].str.extract(r'[ABCDEF](\d{4})') # recupération du numeros de résidu donneur
df_filtre2["id_res_accepteur"] = df_filtre2['Residu2'].str.extract(r'[ABCDEF](\d{4})') # idem accepteur
df_filtre2["chaine"] = df_filtre2["Residu1"].str.extract(r'([ABCDEF])') # récuperation de la chaine
df_filtre2["id_res_doneur"] = df_filtre2["id_res_doneur"].astype(int) #transforme les str des numeros de résidu en int
df_filtre2["id_res_accepteur"] = df_filtre2["id_res_accepteur"].astype(int) #idem


"""Prédiction des structures secondaires grace laisons hydrogène"""

helice = []
beta = []
colonne_result = ["id_res_donneur", "id_res_accepteur", "chaine", "prediction"]
df_result = pd.DataFrame()

for res in range(len(df_filtre2)): #iteration sur l'ensemble des interactions et capture des valeurs
    doneur = df_filtre2["Residu1"][res]
    accepteur = df_filtre2["Residu2"][res]

    doneur_id = df_filtre2["id_doneur"][res]
    accepteur_id = df_filtre2["id_accepteur"][res]
    differrence = abs(doneur_id - accepteur_id)


    if differrence in [3,4,5]: # si interaction de 3 à 5 residu les residu sont ranger dans hélice
        helice.append(doneur_id)
        helice.append(accepteur_id)
        df_result.at[res, "id_res_donneur"] = doneur_id
        df_result.at[res, "id_res_accepteur"] = accepteur_id
        df_result.at[res,"chaine"] = df_filtre2.loc[res,"chaine"]
        df_result.at[res, "prediction"] = "H"
        

    elif 5 > differrence <= 120: #si interaction compris entre 6 et 120 residu les residu sont ranger dans beta
        beta.append(doneur_id)
        beta.append(accepteur)
        df_result.at[res, "id_res_donneur"] = doneur_id
        df_result.at[res, "id_res_accepteur"] = accepteur_id
        df_result.at[res,"chaine"] = df_filtre2.loc[res,"chaine"]
        df_result.at[res, "prediction"] = "B"


"""capture des information des residus de la chaine principales mise dans un dataframe"""

with open('8zap.pdb', 'r') as fichier: #ouverture PDB
    dernier_residu = None
    chaine = None
    residus = []


    for ligne in fichier: # On itere sur toute les ligne puis uniquement aux lignes commençant par 'ATOM'
        
        if ligne.startswith('ATOM'):
            # Extraire les informations pertinentes
            residu_num = ligne[22:26].strip()  # Numéro du résidu (colonnes 23-26)
            residu_nom = ligne[17:20].strip()  # Nom du résidu (colonnes 18-20)
            chaine_id = ligne[21].strip()      # ID de la chaîne (colonne 22)
            
            
            if dernier_residu != residu_num or chaine != chaine_id: # Si on change de résidu ou de chaîne, on fait un traitement
                dernier_residu = residu_num
                chaine = chaine_id
                # Ajouter le nouveau résidu à la liste
                residus.append((residu_num, residu_nom, chaine_id))

df_residus = pd.DataFrame(residus, columns=['Numero_du_residu', 'Nom_du_residu', 'Chaine'])


"""Transfert des données de prédiction du dataframe de hb2 a celui de chaine principale"""

df_result = df_result.reindex(df_residus.index, fill_value= np.nan)
df_residus["prediction"] = ""

for pred in  range(len(df_residus)):

    valeur1 = df_residus.loc[pred, "Numero_du_residu"]
    valeur2 = df_result.loc[pred, "id_res_donneur"]
    valeur3 = df_result.loc[pred, "id_res_accepteur"]
    valeur4 = df_residus.loc[pred, "Chaine"]
    valeur5 = df_result.loc[pred, "chaine"]

    if valeur1 == valeur2 and valeur4 == valeur5:

        df_residus.loc[pred, "prediction"] = df_result.loc[pred, "prediction"]

    if valeur1 == valeur3 and valeur4 == valeur5:
        
        df_residus.loc[pred, "prediction"] = df_result.loc[pred, "prediction"]

"""Lancement du DSSP classique"""

dssp = dssp_originel("8zap.pdb")
dssp_predictions = dssp.read()

"""Comparaison des deux methode"""