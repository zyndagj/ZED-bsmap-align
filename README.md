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

#### Bams
- sample.bam

#### Methylation Calls
- sample_methratio.txt

#### Tiles
- sample_CG.bedgraph
- sample_CHG.bedgraph
- sample_CHH.bedgraph
- sample.tab

#### BigWigs
- sample_CG.bw
- sample_CHG.bw
- sample_CHH.bw

#### YAML
- sample.yaml
