# Meloidogyne spp. genome assemblies

The genome assembly for M. incognita, M. arenaria and M. javanica is described in the pre-print : 

[Unzipped assemblies of polyploid root-knot nematode genomes reveal new kinds of unilateral composite telomeric repeats](https://www.biorxiv.org/content/10.1101/2023.03.29.534350v3.abstract)

The following lines describe the methods used for assign sub-genome for Meloidogyne species, based in synteny and ks values. 

##Separating sub-genomes in Meloidogyne

1- Run McScanX for each species using the default settings

2- Using the collinearity file, run add_kaks_to_MCScanX.pl ([Reubwn](https://github.com/reubwn/collinearity/tree/master))

3- Format the ka/ks file using calculate_collinearity_metric.pl ([Reubwn](https://github.com/reubwn/collinearity/tree/master)) 

4- Select only the synteny groups bigger than 10 genes 

```
python3 groups_bigger_10.py

```
5- Format file "groups_bigger_than_10.txt"

```
awk '$0 !~/#/ {print $3 "\t" $4 "\t" $7}' groups_bigger_than_10.txt > species_values_kaks.txt
```

6- Retrieve the gff file (four columns: contig, gene, start, end) used for McScanX analysis, and run the script on jupyter notebook 

```
Assign_subgenomes_script.py
```

7- To calculate the euclidean distance based on ks values (script by https://github.com/DjampaKozlowski)

```
genome_structure.py
```





