# File `chas_example.txt`
Toy example containing a subset of the alterations of a real ChAS export file.

# File `ascat_example.txt`
Toy example containing a subset of the alterations of a ASCAT result file.

# File `OncoScan_CNV.na33.r2.annot.processed.bed`
Arm coverage of the Oncoscan na33.r2 assay. Of note, the arms 13p, 14p, 15p,
21p and 22p are not covered by the Oncoscan.

## Source
References downloaded from the ThermoFisher website:
OncoScan_CNV_Analysis_Files_NA33.r5.zip. Then the position of each probe was
extracted from the SQLite database (OncoScan_CNV.na33.r2.annot.db file) with
the following command:

`SELECT Chr_id, Start, Stop, Cytoband, dbSNP_RS_ID, ProbeSet_ID FROM
Annotations WHERE Process_Flag > 0;`


Data was then saved in a temporary file (`sql-output.csv`) and processed to 
group probes based on their distance to the next probe (different groups if 
more than 300Kbp, the genome-wide Oncoscan resolution). 
Groups with less than 10 probes and situated at the ends of the
chromosome or close to the centromere were discarded.

Grouping procedure in python3:

```{python}
snp = dict()
with open('sql-output.csv', 'r') as f:
    h = f.readline()  # skip headers
    for l in f:
        toks = l.strip().split('\t')
        arm = toks[0] + toks[3][0]
        pos = int(toks[1])
        if arm not in snp:
            snp[arm] = []
        snp[arm].append(pos)

for arm in sorted(snp.keys()):
    snps = sorted(snp[arm])
    start = {'pos': snps[0], 'i':0}
    for i in range(len(snps)-1):
        if snps[i+1]-snps[i] > 300000:
            print('\t'.join([arm, str(start['pos']), str(snps[i]),
                             str(i+1-start['i']),
                             str((snps[i]-start['pos'])/1000)]))
            start = {'pos': snps[i+1], 'i': i+1}
    print('\t'.join([arm, str(start['pos']), str(snps[-1]),
                     str(len(snps)+1-start['i']),
                     str((snps[-1]-start['pos'])/1000)]))
```

The output of the python script was saved in the file 
`OncoScan_CNV.na33.r2.annot.processed.bed`

## Format
A BED file with 3 columns (chromosome, start, end) and one line for each arm.

# File `LST_gene_list_full_location.txt`
A ChAS export file of a tumor sample with an HR-deficient phenotype

# File `TDplus_gene_list_full_location.txt`
A ChAS export file of a tumor sample with a CDK12-mutated phenotype

# File `triploide_gene_list_full_location.txt`
A ChAS export file of a sample with a triploid tumor
