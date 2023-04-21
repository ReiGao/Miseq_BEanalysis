# Myseq_BEanalysis V1.0

## 1.Raw data cleaning
    Use trimmomatic software to clean the raw data, the parameters are as follows:
```
java -jar trimmomatic-0.36.jar PE -threads $core -phred33  $R1 $R2 $outdir/R1_trim_paired.clean.fq.gz $outdir/R1_trim_unpaired.clean.fq.gz $outdir/R2_trim_paired.clean.fq.gz $outdir/R2_trim_unpaired.clean.fq.gz ILLUMINACLIP:./TruSeq3-PE.fa:2:30:3 LEADING:2 TRAILING:10 SLIDINGWINDOW:1:2 MINLEN:75
```
## 2.Merge PE reads
    Use FLASH software to merge the PE reads, the parameters are as follows:
```
flash -t $core -M 250 -x 0.25 -o $outdir/merge $outdir/R1_trim_paired.clean.fq.gz $outdir/R2_trim_paired.clean.fq.gz
```
## 3.Compute the Baseediting efficiency
    Prepare the barcode with target file named target.txt file,as follows:

```
CGATGT	ACA
CAGATC	ACC
ATCACG	ACG
CATTTT	ACT
CTATAC	CCA
GAGTGG	CCC
GTGAAA	CCG
GTTTCG	CCT
ACTGAT	GCA
CAACTA	GCC
CAGGCG	GCG
ATTCCT	GCT
ACTTGA	TCA
AGTCAA	TCC
ATGTCA	TCG
TAGCTT	TCT
```
    Use the following command to compute the Baseediting efficiency:
```
perl BEanalysisV1.pl ExampleData/example.extendedFrags.fastq ExampleResult> ExampleSample.stat
```



