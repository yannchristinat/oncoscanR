# oncoscanR
Le logiciel OncoscanR est un package R permettant le calcul de différentes caractéristiques propres à une analyse du nombre de copie (CNV) basée sur la technologie Affymetrix Oncoscan. Le logiciel prend en entrée une liste de segments avec leurs nombre de copie mais ne calcule pas les alterations CNV depuis les fichier bruts. OncoscanR suppose que tous les segments qui lui sont donnés sont corrects. Le logiciel calcule la liste des bras chromosomiques qui sont altérés dans leur globalité, le nombre de large-scale state transitions [Popova et al., …], le score TD+ et estime le nombre total de chromosomes (ploïdie). 

Le pipeline standard pour l'Oncoscan se lance avec le Rscript "exec/oncoscan-workflow.R".
Usage:
Rscript path_to_oncoscanR_package/exec/oncoscan-workflow.R CHAS_FILE GENDER
CHAS_FILE: Path to the text export file from ChAS (the software that processes the Oncoscan raw data).
GENDER: Gender of the sample (used to handle sex chromosomes). Has to be M (male) or F (female).

Le script écrit un texte JSON dans le terminal en suivant la structure ci-dessous (exemple):
{
  "armlevel": {
    "AMP": [],
    "LOSS": ["17p", "2q", "4p"],
    "LOH": ["14q", "5q", "8p", "8q"],
    "GAIN": [19p", "19q", "1q", "20p", "20q", "3q", "5p", "6p", "9p", "9q", "Xp", "Xq"]
  },
  "scores": {
    "LST": 12,
    "LOH": 10,
    "TDplus": 66,
    "TD": 18
  },
  "gender": "F",
  "file": "H19001012_gene_list_full_location.txt"
} 
