# projet_court_Gcluzeau M2BI_Projet_court_DSSP

Au cour de ce projet nous avons tenté de réimplémenter le programme de prédiction de structure secondaire DSSP.
Pour ce projet nous debutons d'un fichier PDB à telecharger sur la base de données PDB. Ce fichier est en suite traité avec le programe hbplus.
Celui ci permet de calculer les coordonées des hydrogènes et de calculuer les liaisons hydrogène.

Les liaisons hydrogène de Main chain sont extraite pour prédire les structure secondaires.


## Traiter le format PDB de la protéine par HBPLUS depuis votre terminal. 
### Pour lancer HBPLUS il faut enlever le format .ZIP dans votre terminal.
```bash
unzip hbplus.zip
```

### Se déplacer dans le dossier hbplus
```bash
cd hbplus
```

### Mise en place de l'algorithme 
```bash
make
```

#### Lancer le programme hbplus
```bash
./hbplus
```

## Suivre les instructions du programme hbplus
Le fichier contenant les liaisons hydrogène (fichier .hb2) se trouvera où vous l'avez indiquer.

## Lancer le scipt python


