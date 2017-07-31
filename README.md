# zed-align
Alignment and methylation calling pipeline for Zea Epigenomics Database (ZED).

# Dependencies

* samtools < 1.0.0
  
  From Stampede
  ```shell
  module load samtools/0.1.19
  ```

* BSMAP

  From Stampede
  ```shell
  module load python samtools/0.1.19 bsmap/2.89
  ```

* bedGraphToBigWig

  UCSC Binary Utilities - http://hgdownload.cse.ucsc.edu/admin/exe/

* PyYAML

  http://pyyaml.org/wiki/PyYAML

  From Stampede
  ```
  pip install --user PyYAML
  ```

# Installation

None at the moment. This will change into a proper python module with dependency checking soon enough.

# Usage
```
usage: zed-align.py [-h] -R FASTA -r1 FASTQ [-r2 FASTQ] [-N STR] [-U] [-q Q]
                    [-C Chrom] [-S N] [-d N] [--CG N] [--CHG N] [--CHH N]
                                                                              
Wrapper for Bisulfite Methylation Alignment.                                  
                                                                              
optional arguments:                                                           
  -h, --help          show this help message and exit                         
  -R    FASTA         Reference for alignment                                 
  -r1   FASTQ         Single or first fastq from pair                         
  -r2   FASTQ         Second read                                             
  -O    STR           Output directory (Default: .)
  -N STR, --name STR  Name for run                                            
  -U, --uniq          Only use unique alignments                              
  -q    Q             Fastq Quality Encoding (Default: 33)                    
  -C    Chrom         Chromosome to use for checking bisulfite conversion rate
  -S    N             Window size                                             
  -d    N             Minimum coverage in tile for methylation to be printed  
  --CG  N             Minimum sites per tile (Default: 3)
  --CHG N             Minimum sites per tile (Default: 3)
  --CHH N             Minimum sites per tile (Default: 6)
```

# Output

## Bams
This file fully conforms to BAM/SAM specifications, but has `ZS:Z` and `XR:Z` fields to store strand information. More information can be found at https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation.
### sample.bam

## Methylation Calls
Methylation ratio files generated with BSMAP's [methratio.py](https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation). Field information is listed below, but an in-depth explanation of values can be found at https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation.
### sample_methratio.txt
| Field | Description |
|:-----:|:------------|
| chr | Chromosome |
| pos | 1-index position |
| strand | Chromosome strand |
| context | CG, CHG, or CHH |
| ratio | Effective methylation ratio: C_count / eff_CT_count |
| eff_CT_count | Effective CT coverage: CT*(rev_G/rev_GA) |S
| C_count | Number of C's in reads |
| CT_count | Number of C's or T's in reads |
| rev_G_count | Number of G's on the reverse compliment |
| rev_GA_count | Number of G's or A's on the reverse compliment |
| CI_lower | Lower Wilson score interval of ratio |
| CI_upper | Upper Wilson score interval of ratio |

## Tiles
Tiles of size `S`, with ratios calculated by summing all `C_count` and `CT_count` values from each methylation site in a tile. Whenever a tile is less than the `--CG`, `--CHG`, or `--CHH` criteria, the bedgraph ratios are set to 0, but left unthresholded in the tab file. 0 methylation values are compacted in bedgraph files, but are left expanded in the tab files.

### Example
`sample_methratio.txt`

| chr | pos | strand | context | ratio | eff_CT_count | C_count | CT_count | rev_G_count | rev_GA_count | CI_lower | CI_upper |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 1 | 90 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 94 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 103 | + | CHG | 1.000 | 1.00 | 1 | 1 | 0 | 0 | 0.207 | 1.000 |
| 1 | 113 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 114 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 119 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 125 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 132 | + | CHG | 1.000 | 1.00 | 1 | 1 | 0 | 0 | 0.207 | 1.000 |
| 1 | 133 | + | CG | 1.000 | 1.00 | 1 | 1 | 0 | 0 | 0.207 | 1.000 |
| 1 | 145 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 146 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 149 | + | CG | 1.000 | 1.00 | 1 | 1 | 0 | 0 | 0.207 | 1.000 |
| 1 | 168 | + | CG | 1.000 | 1.00 | 1 | 1 | 0 | 0 | 0.207 | 1.000 |
| 1 | 172 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 177 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 184 | + | CG | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |
| 1 | 187 | + | CHH | 0.000 | 1.00 | 0 | 1 | 0 | 0 | 0.000 | 0.793 |

Looking at the 100-200 bin in this example methratio file, we would have the following values in our tab file.
```
           133 149 168 184
CG_C     = 1 + 1 + 1 + 0
CG_CT    = 1 + 1 + 1 + 1
CG_ratio = 3/4 = 0.75
CG_sites = 4

            103 132
CHG_C     = 1 + 1
CHG_CT    = 1 + 1
CHG_ratio = 2/2 = 1.0
CHG_sites = 2

            113 114 119 125 145 146 172 177 187
CHH_C     = 0 + 0 + 0 + 0 + 0 + 0 + 0 + 0 + 0
CHH_CT    = 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1
CHH_ratio = 0/9 = 0
CHH_sites = 9
```

`sample.tab`

| Chr | Start | End | CG_ratio | CG_C | CG_CT | CG_sites | CHG_ratio | CHG_C | CHG_CT | CHG_sites | CHH_ratio | CHH_C | CHH_CT | CHH_sites |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 1 | 0 | 100 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 2 |
| 1 | 100 | 200 | 0.75 | 3 | 4 | 4 | 1.00 | 2 | 2 | 2 | 0 | 0 | 9 | 9 |

However, since the number of CHG sites is less than 3, the bedgraphs would have the following values.

```
B73_CG.bedgraph.gz
1       0       100     0
1       100     200     0.75
B73_CHG.bedgraph.gz
1       0       1200    0
1       1200    1300    0.56
B73_CHH.bedgraph.gz
1       0       1200    0
1       1200    1300    0.06
```

### Files
- sample_CG.bedgraph
- sample_CHG.bedgraph
- sample_CHH.bedgraph
- sample.tab

## BigWigs
Binary bigwig versions of the bedgraph files.
### Files
- sample_CG.bw
- sample_CHG.bw
- sample_CHH.bw

## YAML
Human and machine readable metadata file about run that also contains simple statistics.
### Files
- sample.yaml
