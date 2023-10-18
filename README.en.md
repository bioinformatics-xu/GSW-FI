# GSW-FI

### Description
{**GSW-FI: A GLM model incorporating shrinkage and double-weighted strategies for identifying cancer driver genes with functional impact**}

### Instructions to GSW-FI
Code description

#### Dependencies
Some R package should be imported to apply GSW-FI, including:

1.  data.table
2.  plyr
3.  stringr
4.  MASS
5.  gamlss

#### Run GSW-FI
Identify driver genes based on GSW-FI

1.  Run data_preprocess.R

- Input: 
    
    maf.file.name

- Output: 
    
    cancer.name + ‘_gene_effect.txt’

    cancer.name + ‘_mutationdata.txt’
   
    cancer.name is the name of the cancer for the to MAF file.

- Folder structure:

    ```
    GSW-FI
    |__ README
    |__ data_preprocess.R
    |__ acc_tcga.maf.txt
    ```

2.  Run calculateFIS_estimatedBFIS.R

- Input: 

    maf.file.name or cancer.name

    cancer.name + ‘_gene_effect.txt’ obtained by data_preprocess.R

    cancer.name + ‘_mutationdata.txt’ obtained by data_preprocess.R

    mutation_type_dictionary_file.txt

    MA_scores_rel3_hg19_full (download from http://mutationassessor.org/r3/)

- Output: 

    cancer.name + ‘_geneFeatureScore.txt’
    
- Folder structure:

    ```
    GSW-FI
    |__ README
    |__ calculateFIS_estimatedBFIS.R
    |__ ACC_gene_effect.txt
    |__ ACC_mutationdata.txt
    |__ mutation_type_dictionary_file.txt
    |__ ./MA_scores_rel3_hg19_full/MA_scores_rel3_hg19_chr
    ```

3.  Run identify_drivers.R

- Input: 

    maf.file.name or cancer.name

    cancer.name + ‘_geneFeatureScore.txt’

- Output: 

    driver genes
    
- Folder structure:

    ```
    GSW-FI
    |__ README
    |__ identify_drivers.R
    |__ ACC_geneFeatureScore.txt
    ```

### Developer

- Xiaolu xu
- lu.xu@mail.dlut.edu.cn 
- School of Computer and Artificial Intelligence
- Liaoning Normal University
- Dalian
- China