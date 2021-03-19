# Projet-Long

Ce logiciel vise à aider à la construction de modèle via I-TASSER de modèle 3D de protéines membranaire. 

## Installation

### HH-blits

La version de hh-blits utilisé est installer par l'environnement conda.  
Téléchargement de la base de donnée PDB70à l'adresse [http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)

### Conda

Ce logiciel à besoins de **python3.5** ou de version plus récentes ainsi que de **conda** ou **miniconda**. 
Merci de vous assurez de ces mises à jours.

### Environnement conda

Une fois ce logiciel cloné, aller dans le répertoire principale
```shell
cd path/to/software
```

Puis lancer la commande :
```shell
conda env create --name envname -f environ.yml
```

## Utilisation

Pensez à activer le shell
```shell
conda activate envname
```
Exemple d'utilisation :
Depuis le répertoire principale du logiciel
```shell
python src/main.py -i path/to/data/fasta.fasta -d path/to/database/pdb70/pdb70 -t 30
```

Depuis n'importe quel repertoire
```shell
python path/to/software/src/main.py -i path/to/data/fasta.fasta -d path/to/database/pdb70/pdb70 -t 30
```
