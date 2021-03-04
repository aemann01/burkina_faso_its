# Raw data processing

## Setup

Your working directory should now be 01-read_processing 

Install the read processing conda environment

```bash
conda env create -f environment.yml
```

To load the environment

```bash
conda activate bf_its_processing
```

Download raw data from ////

```bash
mkdir raw
cd raw
wget ////
cd ..
```

### 1. Install R packages (v3.6.3)

Start an interactive R session and run the following:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")
BiocManager::install("dada2")
BiocManager::install("Biostrings")
BiocManager::install("phyloseq")
# install.packages("tidyverse")
install.packages("stringr")
install.packages("data.table")
install.packages("broom")
install.packages("qualpalr")
install.packages("viridis")
install.packages("seqinr")
install.packages("ape")
install.packages("phytools")
```

```text
Note, if you run into a problem with the magick++ library not installing, try one of the following in a normal terminal:
Configuration failed to find the Magick++ library. Try installing:
 - deb: libmagick++-dev (Debian, Ubuntu)
 - rpm: ImageMagick-c++-devel (Fedora, CentOS, RHEL)
 - csw: imagemagick_dev (Solaris)
 - brew imagemagick@6 (MacOS)
For Ubuntu versions Trusty (14.04) and Xenial (16.04) use our PPA:
   sudo add-apt-repository -y ppa:cran/imagemagick
   sudo apt-get update
   sudo apt-get install -y libmagick++-dev

 If you get the error:

 Cannot find curl-config

 Run (on Ubuntu):

 sudo apt-get install libcurl4-gnutls-dev
 ```

### 2. Load required libraries

```R
library(dada2)
# library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(phyloseq)
library(ape)
library(phytools)
```

### 3. File path setup

```R
rawpath <- "raw/"
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
```

```text
  [1] "EXT-BLK1" "EXT-BLK3" "EXT-BLK4" "EXT-BLK5" "EXT-BLK6" "PCR-BLK2"
  [7] "PCR-BLK3" "PCR-BLK4" "PCR-BLK5" "TM01-01"  "TM01-03"  "TM01-04"
 [13] "TM02-02"  "TM02-03"  "TM02-04"  "TM03-01"  "TM03-02"  "TM03-03"
 [19] "TM03-04"  "TM04-01"  "TM04-02"  "TM04-03"  "TM05-01"  "TM05-02"
 [25] "TM05-03"  "TM06-01"  "TM06-02"  "TM06-03"  "TM07-01"  "TM07-02"
 [31] "TM07-03"  "TM07-04"  "TM08-01"  "TM08-02"  "TM08-03"  "TM09-01"
 [37] "TM09-02"  "TM09-03"  "TM09-04"  "TM10-01"  "TM10-03"  "TM11-02"
 [43] "TM11-04"  "TM12-01"  "TM12-02"  "TM12-03"  "TM12-04"  "TM13-01"
 [49] "TM13-02"  "TM14-01"  "TM14-02"  "TM14-03"  "TM14-04"  "TM15-02"
 [55] "TM16-01"  "TM16-02"  "TM16-03"  "TM16-04"  "TM17-01"  "TM17-02"
 [61] "TM17-03"  "TM17-04"  "TM18-01"  "TM18-02"  "TM18-03"  "TM18-04"
 [67] "TM19-01"  "TM19-02"  "TM19-03"  "TM19-04"  "TM20-01"  "TM20-02"
 [73] "TM20-03"  "TM20-04"  "TM21-01"  "TM21-02"  "TM21-03"  "TM21-04"
 [79] "TM22-03"  "TM22-04"  "TM23-01"  "TM23-02"  "TM23-04"  "TM24-01"
 [85] "TM24-03"  "TM25-01"  "TM25-02"  "TM25-03"  "TM25-04"  "TM26-01"
 [91] "TM26-02"  "TM26-03"  "TM26-04"  "TM27-01"  "TM27-03"  "TM27-04"
 [97] "TM28-01"  "TM28-02"  "TM28-03"  "TM28-04"  "TM29-02"  "TM29-03"
[103] "TM30-01"  "YEAST"
```

### 4. Plot quality scores

```R
system("mkdir img")
png(paste("img/", "forward_quality_plot.png", sep=""))
plotQualityProfile(fnFs[30:35])
dev.off()
png(paste("img/", "reverse_quality_plot.png", sep=""))
plotQualityProfile(fnRs[30:35])
dev.off()
```

![forward quality plot](https://github.com/aemann01/burkina_faso_its/blob/master/01-read_processing/imgs/forward_quality_plot.png)
![reverse quality plot](https://github.com/aemann01/burkina_faso_its/blob/master/01-read_processing/imgs/reverse_quality_plot.png)

The reverse read quality scores are really poor -- just going to run forward reads

### 5. Preliminary filter (removes sequences with N's, truncate to 225bp)

```R
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE, compress = TRUE, truncLen=225)
```

### 6. Primer removal

```R
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
FWD.RC <- dada2:::rc("GCTGCGTTCTTCATCGATGC")
REV.RC <- dada2:::rc("GCTGCGTTCTTCATCGATGC")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "GCTGCGTTCTTCATCGATGC", "-a", REV.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2,"-o", fnFs.cut[i], fnFs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
```

### 7. Quality filter reads

```R
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
out <- filterAndTrim(cutFs, filtFs, minLen=25, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
retained
```

```text
                                   reads.in reads.out percentage_retained
EXT-BLK1_S10_L001_R1_001.fastq.gz       394       385            97.71574
EXT-BLK3_S58_L001_R1_001.fastq.gz        89        72            80.89888
EXT-BLK4_S70_L001_R1_001.fastq.gz        45        43            95.55556
EXT-BLK5_S82_L001_R1_001.fastq.gz       920       870            94.56522
EXT-BLK6_S106_L001_R1_001.fastq.gz      293       258            88.05461
PCR-BLK2_S56_L001_R1_001.fastq.gz       113       104            92.03540
PCR-BLK3_S117_L001_R1_001.fastq.gz      154       141            91.55844
PCR-BLK4_S33_L001_R1_001.fastq.gz       128       121            94.53125
PCR-BLK5_S118_L001_R1_001.fastq.gz      627       563            89.79266
TM01-01_S138_L001_R1_001.fastq.gz       561       481            85.73975
TM01-03_S97_L001_R1_001.fastq.gz       6215      5947            95.68785
TM01-04_S57_L001_R1_001.fastq.gz      10748     10060            93.59881
TM02-02_S139_L001_R1_001.fastq.gz     14648     13938            95.15292
TM02-03_S50_L001_R1_001.fastq.gz       2846      2810            98.73507
TM02-04_S9_L001_R1_001.fastq.gz        4269      4118            96.46287
TM03-01_S123_L001_R1_001.fastq.gz      2670      2588            96.92884
TM03-02_S93_L001_R1_001.fastq.gz        918       879            95.75163
TM03-03_S65_L001_R1_001.fastq.gz        758       610            80.47493
TM03-04_S140_L001_R1_001.fastq.gz      3871      3122            80.65099
TM04-01_S141_L001_R1_001.fastq.gz     17762     17030            95.87884
TM04-02_S116_L001_R1_001.fastq.gz     17617     16944            96.17983
TM04-03_S80_L001_R1_001.fastq.gz      43280     41787            96.55037
TM05-01_S81_L001_R1_001.fastq.gz      21381     19045            89.07441
TM05-02_S27_L001_R1_001.fastq.gz       5604      5533            98.73305
TM05-03_S128_L001_R1_001.fastq.gz      5722      4939            86.31597
TM06-01_S61_L001_R1_001.fastq.gz      10702     10239            95.67371
TM06-02_S89_L001_R1_001.fastq.gz        530       482            90.94340
TM06-03_S109_L001_R1_001.fastq.gz      6679      6568            98.33807
TM07-01_S43_L001_R1_001.fastq.gz       7115      5135            72.17147
TM07-02_S32_L001_R1_001.fastq.gz      16097     15494            96.25396
TM07-03_S105_L001_R1_001.fastq.gz      8086      7139            88.28840
TM07-04_S7_L001_R1_001.fastq.gz        7649      7546            98.65342
TM08-01_S44_L001_R1_001.fastq.gz      63566     60533            95.22858
TM08-02_S31_L001_R1_001.fastq.gz        723       699            96.68050
TM08-03_S101_L001_R1_001.fastq.gz       217       194            89.40092
TM09-01_S104_L001_R1_001.fastq.gz      4136      3907            94.46325
TM09-02_S68_L001_R1_001.fastq.gz      22001     20503            93.19122
TM09-03_S115_L001_R1_001.fastq.gz       271       248            91.51292
TM09-04_S92_L001_R1_001.fastq.gz      30510     30277            99.23632
TM10-01_S42_L001_R1_001.fastq.gz       4031      2223            55.14761
TM10-03_S129_L001_R1_001.fastq.gz      2133      2090            97.98406
TM11-02_S30_L001_R1_001.fastq.gz       2351      2285            97.19268
TM11-04_S78_L001_R1_001.fastq.gz       6364      6086            95.63168
TM12-01_S127_L001_R1_001.fastq.gz      6741      6664            98.85774
TM12-02_S69_L001_R1_001.fastq.gz       5908      5224            88.42248
TM12-03_S67_L001_R1_001.fastq.gz       2579      2085            80.84529
TM12-04_S102_L001_R1_001.fastq.gz      1919      1761            91.76655
TM13-01_S136_L001_R1_001.fastq.gz      8868      8308            93.68516
TM13-02_S6_L001_R1_001.fastq.gz         683       615            90.04392
TM14-01_S41_L001_R1_001.fastq.gz       3536      3275            92.61878
TM14-02_S103_L001_R1_001.fastq.gz      3460      2993            86.50289
TM14-03_S55_L001_R1_001.fastq.gz       9343      8008            85.71123
TM14-04_S98_L001_R1_001.fastq.gz       4972      4914            98.83347
TM15-02_S110_L001_R1_001.fastq.gz      3782      3722            98.41354
TM16-01_S91_L001_R1_001.fastq.gz      12520     12414            99.15335
TM16-02_S85_L001_R1_001.fastq.gz      29885     29444            98.52434
TM16-03_S94_L001_R1_001.fastq.gz        265       265           100.00000
TM16-04_S26_L001_R1_001.fastq.gz       1857      1836            98.86914
TM17-01_S99_L001_R1_001.fastq.gz      16344     16100            98.50710
TM17-02_S112_L001_R1_001.fastq.gz      2703      2619            96.89234
TM17-03_S28_L001_R1_001.fastq.gz       3325      3218            96.78195
TM17-04_S37_L001_R1_001.fastq.gz      52502     49784            94.82305
TM18-01_S39_L001_R1_001.fastq.gz      54281     52864            97.38951
TM18-02_S79_L001_R1_001.fastq.gz      13070     12768            97.68936
TM18-03_S100_L001_R1_001.fastq.gz     41147     40778            99.10322
TM18-04_S8_L001_R1_001.fastq.gz        7017      6843            97.52031
TM19-01_S54_L001_R1_001.fastq.gz        823       803            97.56987
TM19-02_S49_L001_R1_001.fastq.gz      11604     11067            95.37229
TM19-03_S1_L001_R1_001.fastq.gz       10075      9937            98.63027
TM19-04_S124_L001_R1_001.fastq.gz     36775     35864            97.52277
TM20-01_S137_L001_R1_001.fastq.gz     10324      9470            91.72801
TM20-02_S52_L001_R1_001.fastq.gz       1979      1885            95.25013
TM20-03_S63_L001_R1_001.fastq.gz      19534     19209            98.33623
TM20-04_S51_L001_R1_001.fastq.gz        677       639            94.38700
TM21-01_S86_L001_R1_001.fastq.gz       9703      9611            99.05184
TM21-02_S25_L001_R1_001.fastq.gz      15332     15093            98.44117
TM21-03_S77_L001_R1_001.fastq.gz       1513      1284            84.86451
TM21-04_S40_L001_R1_001.fastq.gz      16075     15105            93.96579
TM22-03_S133_L001_R1_001.fastq.gz      5627      5111            90.82993
TM22-04_S64_L001_R1_001.fastq.gz      38749     36803            94.97793
TM23-01_S113_L001_R1_001.fastq.gz       462       434            93.93939
TM23-02_S121_L001_R1_001.fastq.gz     16378     15891            97.02650
TM23-04_S66_L001_R1_001.fastq.gz       6528      5409            82.85846
TM24-01_S76_L001_R1_001.fastq.gz      16759     15505            92.51745
TM24-03_S38_L001_R1_001.fastq.gz       5574      5410            97.05777
TM25-01_S88_L001_R1_001.fastq.gz      78890     78009            98.88326
TM25-02_S111_L001_R1_001.fastq.gz      1929      1838            95.28253
TM25-03_S134_L001_R1_001.fastq.gz      1473      1408            95.58724
TM25-04_S87_L001_R1_001.fastq.gz       8793      8710            99.05607
TM26-01_S74_L001_R1_001.fastq.gz      54689     53289            97.44007
TM26-02_S90_L001_R1_001.fastq.gz       2744      2664            97.08455
TM26-03_S125_L001_R1_001.fastq.gz      3473      3399            97.86928
TM26-04_S114_L001_R1_001.fastq.gz       465       444            95.48387
TM27-01_S126_L001_R1_001.fastq.gz      9093      8840            97.21764
TM27-03_S29_L001_R1_001.fastq.gz       9656      9592            99.33720
TM27-04_S53_L001_R1_001.fastq.gz        735       690            93.87755
TM28-01_S5_L001_R1_001.fastq.gz        4273      3859            90.31126
TM28-02_S4_L001_R1_001.fastq.gz       10843     10587            97.63903
TM28-03_S62_L001_R1_001.fastq.gz       7848      7286            92.83894
TM28-04_S122_L001_R1_001.fastq.gz     12074     11925            98.76594
TM29-02_S135_L001_R1_001.fastq.gz     15465     14997            96.97381
TM29-03_S2_L001_R1_001.fastq.gz       10048      9811            97.64132
TM30-01_S73_L001_R1_001.fastq.gz      57001     54925            96.35796
YEAST_S3_L001_R1_001.fastq.gz          8826      8345            94.55019
```

### 8. Learn and plot error rates

```R
errF <- learnErrors(filtFs, multithread=T, random=T)
png(paste("img/", "error_plot.png", sep=""))
plotErrors(errF, nominalQ=TRUE) 
dev.off()
```

![error plot](https://github.com/aemann01/gut_protozoa/blob/master/01-read_processing/imgs/error_plot.png)

### 9. Dereplication

```R
derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- sample.names
```

### 10. Sample inference

```R
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

### 11. Filter out samples with fewer than 1000 reads post quality filtering

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(derepFs, getN), sapply(dadaFs, getN))
samples_to_keep <- track[,2] > 1000
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)]
paste(samples_to_remove)
```

```text
 [1] "EXT-BLK1" "EXT-BLK3" "EXT-BLK4" "EXT-BLK5" "EXT-BLK6" "PCR-BLK2"
 [7] "PCR-BLK3" "PCR-BLK4" "PCR-BLK5" "TM01-01"  "TM03-02"  "TM03-03" 
[13] "TM06-02"  "TM08-02"  "TM08-03"  "TM09-03"  "TM13-02"  "TM16-03" 
[19] "TM19-01"  "TM20-04"  "TM23-01"  "TM26-04"  "TM27-04" 
```

### 12. Construct sequence table

```R
seqtab <- makeSequenceTable(dadaFs[samples_to_keep])
dim(seqtab)
```

```text
[1]  81 714
```

### 13. Sequence length distribution plot

```R
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))
png(paste("img/", "length_hist.png", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
```

![length histogram](https://github.com/aemann01/gut_protozoa/blob/master/01-read_processing/img/length_hist.png)

### 14. Remove chimeras

```R
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
```

```text
[1]  81 491
```

```R
sum(seqtab.nochim)/sum(seqtab)
```

```text
[1] 0.9141645
```

### 15. Processing summary

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nochimeras")
rownames(track) <- sample.names[samples_to_keep]
track
```

```text
        input filtered denoisedF nochimeras
TM01-03  6215     5947      5923       4991
TM01-04 10748    10060     10042       9993
TM02-02 14648    13938     13899      12453
TM02-03  2846     2810      2789       2686
TM02-04  4269     4118      4109       3904
TM03-01  2670     2588      2565       2417
TM03-04  3871     3122      3102       3076
TM04-01 17762    17030     16983      16582
TM04-02 17617    16944     16912      14699
TM04-03 43280    41787     41764      41203
TM05-01 21381    19045     19027      18667
TM05-02  5604     5533      5514       5434
TM05-03  5722     4939      4919       4910
TM06-01 10702    10239     10198       9977
TM06-03  6679     6568      6510       6469
TM07-01  7115     5135      5091       4833
TM07-02 16097    15494     15461      15326
TM07-03  8086     7139      7132       6700
TM07-04  7649     7546      7497       7487
TM08-01 63566    60533     60483      59814
TM09-01  4136     3907      3890       3208
TM09-02 22001    20503     20469      18934
TM09-04 30510    30277     30260      30120
TM10-01  4031     2223      2206       1906
TM10-03  2133     2090      2078       1951
TM11-02  2351     2285      2269       2228
TM11-04  6364     6086      6060       5819
TM12-01  6741     6664      6652       5930
TM12-02  5908     5224      5183       4969
TM12-03  2579     2085      2070       1659
TM12-04  1919     1761      1747       1703
TM13-01  8868     8308      8280       7891
TM14-01  3536     3275      3261       2104
TM14-02  3460     2993      2983       2765
TM14-03  9343     8008      7979       7979
TM14-04  4972     4914      4904       4839
TM15-02  3782     3722      3698       3677
TM16-01 12520    12414     12400      11474
TM16-02 29885    29444     29411      27265
TM16-04  1857     1836      1821       1451
TM17-01 16344    16100     16080      13935
TM17-02  2703     2619      2584       2385
TM17-03  3325     3218      3202       2018
TM17-04 52502    49784     49717      44146
TM18-01 54281    52864     52798      48932
TM18-02 13070    12768     12705      11541
TM18-03 41147    40778     40689      39143
TM18-04  7017     6843      6828       6354
TM19-02 11604    11067     11029      10958
TM19-03 10075     9937      9924       7308
TM19-04 36775    35864     35835      34537
TM20-01 10324     9470      9437       7986
TM20-02  1979     1885      1868       1758
TM20-03 19534    19209     19172      18719
TM21-01  9703     9611      9590       7048
TM21-02 15332    15093     15063      11644
TM21-03  1513     1284      1278        588
TM21-04 16075    15105     15064      10480
TM22-03  5627     5111      5076       4489
TM22-04 38749    36803     36768      32287
TM23-02 16378    15891     15824      13193
TM23-04  6528     5409      5368       3626
TM24-01 16759    15505     15470       9731
TM24-03  5574     5410      5376       4946
TM25-01 78890    78009     77982      72167
TM25-02  1929     1838      1791       1697
TM25-03  1473     1408      1379       1322
TM25-04  8793     8710      8696       8402
TM26-01 54689    53289     53268      46880
TM26-02  2744     2664      2653       2156
TM26-03  3473     3399      3387       3105
TM27-01  9093     8840      8835       5661
TM27-03  9656     9592      9575       8848
TM28-01  4273     3859      3834       3419
TM28-02 10843    10587     10562       9464
TM28-03  7848     7286      7249       5584
TM28-04 12074    11925     11914      11831
TM29-02 15465    14997     14951      13781
TM29-03 10048     9811      9794       9439
TM30-01 57001    54925     54884      49387
YEAST    8826     8345      8315       8260
```

### 16. Save output

```R
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
```

### 17. Clean up ASV names

```R
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")
```

### 18. Assign taxonomy using BLAST and MEGAN6

Download nt blast database

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*tar.gz
ls *gz | parallel 'tar xzf {}'
```

Run sequences through BLAST

```bash
blastn -db /home/allie/refdb/nt/nt -query rep_set.fa -out rep_set.blast.out -evalue 1e-10 -outfmt 6
```

Get taxonomic assignments with MEGAN6 and clean up

```bash
sed 's/;/\t/g' rep_set.tax.txt | sed 's/"//g' > temp
mv temp rep_set.tax.txt
```

Next need to merge ASV and taxonomy information (ASVs with no hits are dropped).

```R
seq.tab <- read.table("sequence_table.merged.txt", sep="\t", header=T, row.names=1)
tax.tab <- read.table("rep_set.lowest-tax.txt", sep="\t", header=T, row.names=1)
tran.seq <- t(seq.tab)
merge.tab <- merge(tran.seq, tax.tab, by=0, all=TRUE)
write.table(merge.tab, file="sequence_table.taxonomy.txt", quote=F, sep="\t")
```



## Metagenomic data processing

### 1. Pull data (from Jacobson et al. 2021)

First need to get data that matches our samples (n=59)

```bash
mkdir raw_WGS
cd raw_WGS
prefetch --option-file ../wgs.query -O .
mv SRR13378*/*sra .
find -empty -type d -delete
ls *sra | parallel 'fasterq-dump --split-files {}'
rm *sra
```

### 2. Trim adapters, quality filter, and merge paired end reads

First need to identify the adapter sequence

```bash
AdapterRemoval --identify-adapters --file1 SRR13378546_1.fastq --file2 SRR13378546_2.fastq
```

Trim off adapters

```bash
ls *_1* | sed 's/_1.fastq//' | parallel 'AdapterRemoval --file1 {}_1.fastq --file2 {}_2.fastq --trimns --trimqualities --minquality 30 --gzip --collapse --basename {} --minlength 100 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
```

Trim poly G tails

cutadapt -a "A{100}" -o output.fastq input.fastq

### 3. Dereplication

Keep collapsed, singletons, high quality forward reads

```bash
rm *pair2* *singleton* *discarded*
ls *settings | sed 's/.settings//' | parallel 'cat {}.collapsed.gz {}.pair1.truncated.gz {}.collapsed.truncated.gz > {}.full.gz'
rm *collapsed* *pair1* *fastq
```

Get a file of accession numbers and original sample IDs

```bash
awk -F"\t" '{print $1, "\t", $33}' sra.info | sed 's/ //g' > wgs.rename
```

Convert to fasta

```bash
ls *full.gz | parallel 'gzip -d {}'
ls *full | parallel 'seqtk seq -a {} > {}.fa'
rm *full
```

Dereplicate

```bash
ls *fa | sed 's/.full.fa//' | while read line; do vsearch --derep_fulllength $line.full.fa --output $line.uniq.fa; done
rm *full.fa
```



Change sequence headers to sample IDs

```bash
tagSeq()
{
    set $file
    for i in $name; do
        sed "s/^>/>${i}_/" ${1} > ${1}_fix
        shift
    done
}
name=$(awk '{print $2}' ../wgs.rename)
file=$(awk '{print $1}' ../wgs.rename)
tagSeq
rm *fna
cat *fix > mock_oral.fa
rm *fix
```









Concatenate all samples and clean up

```bash
cat *fa > all_samples.fa
rm *full* *gz 
```

### 4. Assign taxonomy with KRAKEN2

Build and index database (only run once, takes a long time to complete)

```R
system("kraken2-build --standard --db nt")
```

Assign taxonomy using KRAKEN2

```R
kraken2 --db ~/kraken_nt/ --threads 8 --use-names --output all_samples.uniq.tax all_samples.uniq.fa.gz
```









