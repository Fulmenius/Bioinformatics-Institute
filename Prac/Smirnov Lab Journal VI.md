## 0.Obtaining the data
```shell
wget -o ./raw_data/ferm_0_min_rep_1.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941816/SRR941816.fastq.gz  
... 
And so on
```

## 1. Aligning with HISAT2
The index was built with
```shell
hisat2-build reference_genome.fna genome_index.sam 
```
The alignment was performed with the following commands:
```shell
hisat2 -p 4 -x genome_index.sam -U ./ferm_0_min_rep_1.fastq | samtools sort > out_time_zero_1.bam

hisat2 -p 4 -x genome_index.sam -U ./ferm_0_min_rep_2.fastq | samtools sort > out_time_zero_2.bam

hisat2 -p 4 -x genome_index.sam -U ./ferm_30_min_rep_1.fastq | samtools sort > out_time_30_1.bam

hisat2 -p 4 -x genome_index.sam -U ./ferm_30_min_rep_2.fastq | samtools sort > out_time_30_2.bam
```

## 2. Determining the expression levels
With the help of `gffread` the [[.gff]] annotation file was converted to [[.gtf]], after which [[featureCounts]] was used to determine expression levels. The resulting files were simplified by discarding columns 1 and 7-10. 
```shell
featureCounts -g gene_id -a annotation_file.gtf -o feature_counts out_time_*
cat feature_counts | cut -f 1,7-10 > simple_counts.txt 
```

## 3. Finding the difference
The features differing between the 0 min and the 30 min sequences were found with the help of [[Deseq2]]. 
```shell
cat simple_counts.txt | R -f deseq2.r 
cat norm-matrix-deseq2.txt | R -f draw-heatmap.r  
```

## 4. Interpreting results
We use 
```shell
head -n 50 result.txt | cut -f 1 | cut -d "-" -f 2 > genes.txt
```
to leave only Top-50 genes sorted by [[Adjusted p-value|adjusted p-value]]. 
We then thried to use materials available by the link: http://www.yeastgenome.org/cgi-bin/GO/goSlimMapper.pl to analyze the genes' functions, but the app didn't work, so we resorted to https://go.princeton.edu/cgi-bin/GOTermMapper instead. The results were the following (names of individual genes omitted, only categories containing  a non-zero number of genes are shown)?
| GO Term                                                                                               | GO Term Usage in Gene List | Genome Frequency of Use             |
| ----------------------------------------------------------------------------------------------------- | -------------------------- | ----------------------------------- |
| [ribosome biogenesis](http://amigo.geneontology.org/amigo/term/GO:0042254)                            | 20 of 47 genes, 42.55%     | 491 of 6464 annotated genes, 7.60%  |
| [protein-containing complex assembly](http://amigo.geneontology.org/amigo/term/GO:0065003)            | 10 of 47 genes, 21.28%     | 574 of 6464 annotated genes, 8.88%  |
| [transmembrane transport](http://amigo.geneontology.org/amigo/term/GO:0055085)                        | 7 of 47 genes, 14.89%      | 491 of 6464 annotated genes, 7.60%  |
| [mRNA metabolic process](http://amigo.geneontology.org/amigo/term/GO:0016071)                         | 7 of 47 genes, 14.89%      | 371 of 6464 annotated genes, 5.74%  |
| [DNA-templated transcription](http://amigo.geneontology.org/amigo/term/GO:0006351)                    | 5 of 47 genes, 10.64%      | 286 of 6464 annotated genes, 4.42%  |
| [carbohydrate metabolic process](http://amigo.geneontology.org/amigo/term/GO:0005975)                 | 4 of 47 genes, 8.51%       | 291 of 6464 annotated genes, 4.50%  |
| [amino acid metabolic process](http://amigo.geneontology.org/amigo/term/GO:0006520)                   | 4 of 47 genes, 8.51%       | 246 of 6464 annotated genes, 3.81%  |
| [regulation of DNA-templated transcription](http://amigo.geneontology.org/amigo/term/GO:0006355)      | 3 of 47 genes, 6.38%       | 594 of 6464 annotated genes, 9.19%  |
| [tRNA metabolic process](http://amigo.geneontology.org/amigo/term/GO:0006399)                         | 3 of 47 genes, 6.38%       | 221 of 6464 annotated genes, 3.42%  |
| [carbohydrate derivative metabolic process](http://amigo.geneontology.org/amigo/term/GO:1901135)      | 3 of 47 genes, 6.38%       | 360 of 6464 annotated genes, 5.57%  |
| [nucleocytoplasmic transport](http://amigo.geneontology.org/amigo/term/GO:0006913)                    | 3 of 47 genes, 6.38%       | 191 of 6464 annotated genes, 2.95%  |
| [protein catabolic process](http://amigo.geneontology.org/amigo/term/GO:0030163)                      | 2 of 47 genes, 4.26%       | 316 of 6464 annotated genes, 4.89%  |
| [intracellular protein transport](http://amigo.geneontology.org/amigo/term/GO:0006886)                | 2 of 47 genes, 4.26%       | 460 of 6464 annotated genes, 7.12%  |
| [lipid metabolic process](http://amigo.geneontology.org/amigo/term/GO:0006629)                        | 2 of 47 genes, 4.26%       | 361 of 6464 annotated genes, 5.58%  |
| [sulfur compound metabolic process](http://amigo.geneontology.org/amigo/term/GO:0006790)              | 1 of 47 genes, 2.13%       | 134 of 6464 annotated genes, 2.07%  |
| [DNA replication](http://amigo.geneontology.org/amigo/term/GO:0006260)                                | 1 of 47 genes, 2.13%       | 165 of 6464 annotated genes, 2.55%  |
| [signaling](http://amigo.geneontology.org/amigo/term/GO:0023052)                                      | 1 of 47 genes, 2.13%       | 447 of 6464 annotated genes, 6.92%  |
| [cytoplasmic translation](http://amigo.geneontology.org/amigo/term/GO:0002181)                        | 1 of 47 genes, 2.13%       | 204 of 6464 annotated genes, 3.16%  |
| [generation of precursor metabolites and energy](http://amigo.geneontology.org/amigo/term/GO:0006091) | 1 of 47 genes, 2.13%       | 227 of 6464 annotated genes, 3.51%  |
| [snRNA metabolic process](http://amigo.geneontology.org/amigo/term/GO:0016073)                        | 1 of 47 genes, 2.13%       | 34 of 6464 annotated genes, 0.53%   |
| [protein modification process](http://amigo.geneontology.org/amigo/term/GO:0036211)                   | 1 of 47 genes, 2.13%       | 839 of 6464 annotated genes, 12.98% |
| [DNA recombination](http://amigo.geneontology.org/amigo/term/GO:0006310)                              | 1 of 47 genes, 2.13%       | 268 of 6464 annotated genes, 4.15%  |
The table is available at https://go.princeton.edu/tmp//18186_slimTerms.html. 
