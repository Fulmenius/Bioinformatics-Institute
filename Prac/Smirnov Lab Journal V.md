### Tasks:
* Establish maternal ([[DNA#Mytochodrial DNA|mtDNA]]) and paternal ([[Y-chromosome|Y-chromosome]]) [[Haplogroups|haplogroups]], and, optionally, probable ethnicity. 
* Annotate the obtained [[SNP - Single Nucleotide Polymorphisms|SNPs]] and extract all clinically relevant SNPs from the [[ClinVar]] database. Make a [[.vcf|.vcf]] file from the raw reads. Identify phenotypic traits from the raw data ("genome sketching")
* Make changes to a given genome to get a person with the desired characteristics. I.e., find 5 (well, yes, but actually 10) variants that would be useful, and specify what nucleotide exactly and in what position in this genome needs to be replaced using [[CRISPR]]. 

## Data obtention
The data was obtained by the [link](https://drive.google.com/file/d/1QJkwJe5Xl_jSVpqdTSNXP7sqlYfI666j/view) and unpacked using `unzip`. 

## File conversion
[[PLINK]] v1.9 was used for file conversion to the .vcf format and initial processing.
All SNPs corresponding to deletions and insertions were removed to make the file compatible with annotation tools:
```shell
plink1.9 --23file ./SNP_raw.txt --recode vcf --out snps_clean --output-chr MT --snps-only just-acgt
```

## Establishing haplogroups
The paternal and maternal haplogroups were established.

### mtDNA Haplogroup
mtDNA haplogroup was established using the online tool https://dna.jameslick.com/mthap/mthap.cgi. mthap version 0.19b (2015-05-11)
The best mtDNA haplogroup match was H(T152C) of the [H haplogroup](https://ru.wikipedia.org/wiki/%D0%93%D0%B0%D0%BF%D0%BB%D0%BE%D0%B3%D1%80%D1%83%D0%BF%D0%BF%D0%B0_H_(%D0%BC%D1%82%D0%94%D0%9D%D0%9A)). 

### Y-DNA Haplogroup
https://ytree.morleydna.com was used to predict the most likely Y-chromosome subclade. The result was [R1a1a](https://ru.wikipedia.org/wiki/%D0%93%D0%B0%D0%BF%D0%BB%D0%BE%D0%B3%D1%80%D1%83%D0%BF%D0%BF%D0%B0_R1a1_(Y-%D0%94%D0%9D%D0%9A)). 

Based on this results one can hardly conclude anything specific about the ethnicity of the professor. 

### For future work: PCA analysis
PCA analysis was not performed both because of the lack of time and the lack of disk space. 
For future reference the instruction is [here](https://www.biostars.org/p/335605/). 

## Annotation: **♂**sex**♂** and eye color

## Sex identification
Y-chromosome is present, so the most probable variant is male. 

## Eye color identification
![](https://lh6.googleusercontent.com/5IZDYSxzg77DxJfEp1CQTRTnpYHAKfNxCZEp6GV5b50-d4xR1tZ1DB8TyBwfoIYqqZdDtJ5JnwN1iXLheUWBEkXSjOsPSzKZPQCr-YmdCY-_B3Y6uZPqcr3A3fxebxa2zrCf3ALVf9AJDzXT2_wDycA)
![](https://lh5.googleusercontent.com/98yiuunmly_iQYeXo8uGnInHu_w8Vke_fYiKmmHrHPcJi4ljGAkln2-q9mCCFfDwlalDEFDZI0C5xSnIAyUO8QzsuxTQTj0xCtiNac35grQrZh7cPtw7hbj0pdYEXpUKxbN93C8xbDHdhNtXZDOV0sE)
The eye and the skin colors were estimated based on the algorithm proposed in [this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694299/) article (see the illustration above). To analyze the SNPs, [[IGV browser|IGV browser]] was used. The .vcf file was also assessed manually to determine if 1/1 SNPs were present. The following SNPs are observed:
| Position   | SNP                        | Cocnlusion                                                                                        |
| ---------- | -------------------------- | ------------------------------------------------------------------------------------------------- |
| rs12913832 | Reference: A, Alternate: G | A/G Not Blue                                                                                      |
| rs16891982 | Reference: C, Alternate: G | C/G No conclusion                                                                                 |
| rs6119471  | Not Present                | ---                                                                                               |
| rs12203592 | Reference: C, Alternate: T | No conclusion                                                                                     |
| rs12896399 | Reference: G, NO VARIATION | G/G No conlusion based on the diagram, but if rs12913832 were A/A, this would indicate brown eyes |
|            |                            |                                                                                                   |

## Annotation of all SNPs, selection of clinically relevant
```shell
Was performed by my collegue because of a software problem
```

## Improvements
### Part I. Fixing what may harm

### Part II. TRANSHUMANISM!!!
| Position                                                 | Goal                            | Current SNP or MNP/ Proposed modification                                                 | Possible result                                                                                   |
| -------------------------------------------------------- | ------------------------------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| [rs5882](https://www.snpedia.com/index.php/Rs5882)       | Longevity, Health, Intelligence | AA/GG                                                                                       | Lower risk of Alzheimer's disease and dementia, higher chances for longevity                      |
| [rs363050](https://www.snpedia.com/index.php/Rs363050) |  Intelligence |     AG/AA                      | +2.84 PIQ on average |
| [rs2811712](https://www.snpedia.com/index.php/Rs2811712) | Health, Longevity  |AA/GG  |    1.5x less risk for physical impairment with age.    |
| [rs2802292](https://www.snpedia.com/index.php/Rs2802292) | Longevity                       | GT / GG                                                                                      | Higher chance of extreme longevity                                                                |
| [rs16964465](https://www.snpedia.com/index.php/Rs16964465)      |  Health             |    AA/CC      |       Possible resistance to obesity                                                                    |
