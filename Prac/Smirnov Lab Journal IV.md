## 1. Obtaining data
The files were downloaded with the use of links in the [assignment text](https://docs.google.com/document/d/135503Q9jSzbNBHl_fQZXK21l3O4ZUDMV8OuOPLCr-6U/edit):
```shell
wget -O getAnnoFasta.pl http://augustus.gobics.de/binaries/scripts/getAnnoFasta.pl
```
etc. 

## 2. Structural annotation
The [[Proteins|protein]] sequences were obtained with 
```shell
perl getAnnoFasta.pl augustus.whole.gff
```
The number of proteis was determined with the use of 
```shell
grep '>' augustus.whole.aa | wc -l
```
as 16435.

## 3. Physical localization
The list of peptides obtained through [[Tandem mass spectroscopy|tandem mass spectrometry]] was downloaded manually through graphic interface because of the author's lack of understanding of both Yandex API and Linux. 
The standard [[BLAST]]+ package was used to align the tandem mass spectrometry results on the  peptide reference:
```shell
makeblastdb -in augustus.whole.aa -dbtype prot -out blast_database
```

```shell
blastp -db blast_database -query peptides.fa -outfmt 6 -out blast_output_alignment
```

The result contains 107 positions. If we only take the unique names, we get a list of 34 positions:
```shell
awk '{print $2}' blast_output_alignment | sort | uniq > search_in_aa
```

```shell
g10513.t1
g10514.t1
g11320.t1
g11513.t1
g11806.t1
g11960.t1
g12388.t1
g12510.t1
g12562.t1
g1285.t1
g13530.t1
g14472.t1
g15153.t1
g15484.t1
g16318.t1
g16368.t1
g2203.t1
g3428.t1
g3679.t1
g4106.t1
g4970.t1
g5237.t1
g5443.t1
g5467.t1
g5502.t1
g5503.t1
g5510.t1
g5616.t1
g5641.t1
g5927.t1
g702.t1
g7861.t1
g8100.t1
g8312.t1
```

Using _[[samtools]] faidx_, we extract these sequences from the original file in [[FASTA]] format:
```shell
for f in (cat search_in_aa); samtools faidx augustus.whole.aa $f >> test_cycle_augustus; end
```

## 4. Localization prediction

### 4a. WoLF PSORT
We use the `test_cycle_augustus` file obtined at the previous step as an input for the [WoLF PSORT tool](https://wolfpsort.hgc.jp/); we obtain
```css
g10513.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed1.html#g10513.t1) nucl: 20, cyto_nucl: 14.5, cyto: 7, extr: 3, E.R.: 1, golg: 1 g10514.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed2.html#g10514.t1) nucl: 19, cyto_nucl: 15, cyto: 9, extr: 3, mito: 1 g11320.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed3.html#g11320.t1) plas: 24.5, extr_plas: 16, extr: 6.5, lyso: 1 g11513.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed4.html#g11513.t1) cyto: 17, cyto_nucl: 12.8333, cyto_mito: 9.83333, nucl: 7.5, E.R.: 3, mito: 1.5, plas: 1, pero: 1, golg: 1 g11806.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed5.html#g11806.t1) nucl: 18, cyto_nucl: 11.8333, mito: 5, extr: 4, cyto: 3.5, cyto_pero: 2.66667, cysk_plas: 1 g11960.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed6.html#g11960.t1) nucl: 32 g12388.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed7.html#g12388.t1) extr: 25, plas: 4, mito: 2, lyso: 1 g12510.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed8.html#g12510.t1) plas: 29, cyto: 3 g12562.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed9.html#g12562.t1) extr: 30, lyso: 2 g1285.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed10.html#g1285.t1) extr: 25, plas: 5, mito: 1, lyso: 1 g13530.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed11.html#g13530.t1) extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5 g14472.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed12.html#g14472.t1) nucl: 28, plas: 2, cyto: 1, cysk: 1 g15153.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed13.html#g15153.t1) extr: 32 g15484.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed14.html#g15484.t1) nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1 g16318.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed15.html#g16318.t1) nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1 g16368.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed16.html#g16368.t1) nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1 g2203.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed17.html#g2203.t1) plas: 29, nucl: 2, golg: 1 g3428.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed18.html#g3428.t1) mito: 18, cyto: 11, extr: 2, nucl: 1 g3679.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed19.html#g3679.t1) extr: 26, mito: 2, lyso: 2, plas: 1, E.R.: 1 g4106.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed20.html#g4106.t1) E.R.: 14.5, E.R._golg: 9.5, extr: 7, golg: 3.5, lyso: 3, pero: 2, plas: 1, mito: 1 g4970.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed21.html#g4970.t1) plas: 32 g5237.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed22.html#g5237.t1) plas: 24, mito: 8 g5443.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed23.html#g5443.t1) extr: 28, nucl: 3, cyto: 1 g5467.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed24.html#g5467.t1) extr: 27, plas: 4, mito: 1 g5502.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed25.html#g5502.t1) extr: 31, lyso: 1 g5503.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed26.html#g5503.t1) extr: 29, plas: 1, mito: 1, lyso: 1 g5510.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed27.html#g5510.t1) plas: 23, mito: 7, E.R.: 1, golg: 1 g5616.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed28.html#g5616.t1) extr: 31, mito: 1 g5641.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed29.html#g5641.t1) extr: 31, lyso: 1 g5927.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed30.html#g5927.t1) nucl: 30.5, cyto_nucl: 16.5, cyto: 1.5 g702.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed31.html#g702.t1) extr: 29, plas: 2, lyso: 1 g7861.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed32.html#g7861.t1) nucl: 16, cyto_nucl: 14, cyto: 8, plas: 5, pero: 1, cysk: 1, golg: 1 g8100.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed33.html#g8100.t1) nucl: 16.5, cyto_nucl: 12.5, cyto: 7.5, plas: 5, extr: 2, E.R.: 1 g8312.t1 [details](https://wolfpsort.hgc.jp/results/aST735b0c4fbde0e36c16a564333ebef272.detailed34.html#g8312.t1) nucl: 15.5, cyto_nucl: 15.5, cyto: 12.5, mito: 2, plas: 1, golg: 1
```

### 4b. TargetP 1.1 Server
The same file was given to the [TargetP 2.0 online tool](https://services.healthtech.dtu.dk/service.php?TargetP-2.0). The results are:
```python
# TargetP-2.0	Organism: Non-Plant	Timestamp: 20221215230923
# ID	Prediction	OTHER	SP	mTP	CS Position
g10513.t1	OTHER	0.999999	0.000001	0.000000	
g10514.t1	OTHER	0.999543	0.000349	0.000107	
g11320.t1	SP	0.000184	0.999816	0.000000	CS pos: 20-21. AYS-AG. Pr: 0.7236
g11513.t1	OTHER	0.999434	0.000401	0.000164	
g11806.t1	OTHER	0.998977	0.000887	0.000136	
g11960.t1	OTHER	0.999996	0.000002	0.000002	
g12388.t1	SP	0.000490	0.999481	0.000029	CS pos: 16-17. ASA-SS. Pr: 0.6485
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g1285.t1	SP	0.003029	0.996798	0.000173	CS pos: 16-17. ASA-TS. Pr: 0.7127
g13530.t1	SP	0.116007	0.883840	0.000153	CS pos: 19-20. TIP-FT. Pr: 0.3552
g14472.t1	OTHER	0.999999	0.000001	0.000000	
g15153.t1	SP	0.000014	0.999986	0.000000	CS pos: 16-17. AYA-AN. Pr: 0.8378
g15484.t1	OTHER	0.999980	0.000010	0.000010	
g16318.t1	OTHER	0.997047	0.002953	0.000000	
g16368.t1	OTHER	0.996693	0.003307	0.000000	
g2203.t1	OTHER	0.999869	0.000031	0.000100	
g3428.t1	OTHER	0.999903	0.000033	0.000064	
g3679.t1	SP	0.001755	0.998023	0.000222	CS pos: 18-19. TFA-AR. Pr: 0.5523
g4106.t1	OTHER	0.729658	0.266917	0.003425	
g4970.t1	OTHER	0.999996	0.000003	0.000001	
g5237.t1	OTHER	0.999545	0.000345	0.000111	
g5443.t1	OTHER	0.952853	0.043784	0.003363	
g5467.t1	SP	0.000096	0.999845	0.000059	CS pos: 16-17. ASA-GS. Pr: 0.6543
g5502.t1	SP	0.001134	0.998823	0.000043	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5503.t1	SP	0.001222	0.998720	0.000058	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5510.t1	OTHER	0.999108	0.000016	0.000876	
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g5927.t1	OTHER	0.999995	0.000001	0.000004	
g702.t1	SP	0.000347	0.999652	0.000001	CS pos: 16-17. ALA-AN. Pr: 0.8153
g7861.t1	OTHER	0.999975	0.000004	0.000022	
g8100.t1	OTHER	0.999955	0.000024	0.000021	
g8312.t1	OTHER	0.999930	0.000065	0.000004	
```

## 5. BLAST search
The [[BLAST]] search was performed by my collaborator. The results are included in the final table.

## 6. HMM search
HMM search results:
| Peptide name                                                        | HMM hits                                                                                                                   | E-values (Ind) |
| ------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | -------------- |
| g11513.t1                                                           | Transport protein Trs120 or TRAPPC9, TRAPP II complex subunit                                                              | 9.6e-10        |
| g11960.t1 / [zf-C3HC4](http://pfam.xfam.org/family/zf-C3HC4)        | Zinc finger, C3HC4 type (RING finger)                                                                                      | 4.2e-05        |
| g12388.t1 / CBM_14                                                  | Chitin binding Peritrophin-A domain  (multiple sequences of this type are present, only this is aknowledged in this table) | 2.5e-06        |
| g15484.t1 / Vps51                                                   | Vps51/Vps67                                                                                                                |                |
| g2203.t1 / Glyco_hydro_31                                           | Glycosyl hydrolases family 31                                                                                              | 4.8e-41        |
| g3428.t1 / EF-hand_1                                                | EF hand                                                                                                                    | 4.8e-10        |
| g3679.t1 / Astacin                                                  | Astacin (Peptidase family M12A)                                                                                            | 2.6e-09        |
| g4970.t1 / Trypsin                                                  | Trypsin                                                                                                                    | 1.9e-21        |
| g5510.t1 / MARVEL                                                   | Membrane-associating domain (SERIOUSLY? EVEN IN THE TARDIGRADES GENOME? THE STUDIO BOSSES MUST BE NUTS)                    | 1.8e-09 \      |
| g7861.t1 / [SNF2-rel_dom](http://pfam.xfam.org/family/SNF2-rel_dom) | SNF2-related domain                                                                                                        | 1.2e-28        |
| g8100.t1 / [Inositol_P](http://pfam.xfam.org/family/Inositol_P)     | Inositol monophosphatase family                                                                                            | 1.9e-37        |
| g8312.t1 / [Clathrin](http://pfam.xfam.org/family/Clathrin)         | Region in Clathrin and VPS                                                                                                 | 5.4e-23        |

## 7. Integration and differentiation
The table can be found by the following address: https://docs.google.com/spreadsheets/d/1DonINJTRCSwUB_vSRdA0TEN1JWgnoNqn2iupZ-6EbJE/edit#gid=0
