# Rainbow Trout 0dph RRBS Analysis via methylKit"
# Updated: 15.12.2021

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Setup
```{r, messages=FALSE, results=FALSE}
# setwd("/data/user/kfreij95/rrbs_trout_kwf/rrbs_files")

# browseVignettes("methylKit")
# package.version("methylKit")

# Github page for methylKit: https://github.com/al2na/methylKit

library("methylKit")
library("vegan")
library("gclus")
```
## Key

1. Treatment A(1) - Choline Deficient Diet (Chol-) - 11 replicates
2. Treatment B(0) - Choline Adequate Diet (Chol) - 10 replicates
3. Treatment C(2) - Choline Supplemented Diet (Chol+) - 14 replicates

# Reading in Data
```{r reading in data, message=FALSE, echo=FALSE, eval=FALSE, results=FALSE}
#Sort your SAM files
##via sort_sam.sh in sam_files/

# Reading one bismark file:
#my.file=system.file("extdata", "test.fastq_bismark.sorted.min.sam",
#                                                       package = "methylKit")
#obj=processBismarkAln(my.file,"test",assembly="hg18",save.folder=NULL,
#                 save.context="CpG",read.context="CpG")

# Reading multiple files
file.list2=list("/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/32_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/40_S2_R1_001_bismark_bt2_pe.sort.sam", 
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/41_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/43_S16_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/46_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/47_S5_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/50_S6_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/54_S7_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/59_S8_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/62_S9_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/64_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/76_S2_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/77_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/78_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/79_S5_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/81_S6_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/82_S7_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/85_S8_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/86_S9_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/90_S10_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/95_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/96_S2_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/99_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/100_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/101_S5_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/102_S6_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/108_S7_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/109_S8_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/117_S9_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/118_S10_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/120_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/121_S2_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/124_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/125_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/127_S5_R1_001_bismark_bt2_pe.sort.sam")

file.list3=list()

# Double check the files exist in the working directory
file.exists("/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/32_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/40_S2_R1_001_bismark_bt2_pe.sort.sam", 
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/41_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/43_S16_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/46_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/47_S5_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/50_S6_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/54_S7_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/59_S8_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/62_S9_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/64_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/76_S2_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/77_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/78_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/79_S5_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/81_S6_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/82_S7_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/85_S8_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/86_S9_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/90_S10_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/95_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/96_S2_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/99_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/100_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/101_S5_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/102_S6_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/108_S7_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/109_S8_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/117_S9_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/118_S10_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/120_S1_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/121_S2_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/124_S3_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/125_S4_R1_001_bismark_bt2_pe.sort.sam",
            "/data/user/kfreij95/rrbs_trout_kwf/rrbs_files/sam_files/127_S5_R1_001_bismark_bt2_pe.sort.sam")

# Read Bismark files (sorted SAM files) --> methylRaw/methylRawList
# Can have more than one treatment vector as seen below.
"0dph_rrbsDB"=processBismarkAln(file.list2,
            sample.id=list("32_C","40_C","41_A","43_A","46_B","47_C","50_B","54_A","59_C",
                           "62_A","64_C","76_B", "77_B", "78_C","79_B","81_C","82_A",
                           "85_B","86_C","90_B","95_B","96_C","99_C", "100_C","101_A",
                           "102_A","108_A","109_B","117_A","118_A","120_B","121_C","124_C",
                           "125_C","127_A"), 
            assembly="omyk1.1A",
            treatment=c(2,2,1,1,0,2,0,1,2,
                        1,2,0,0,2,0,2,1,
                        0,2,0,0,2,2,2,1,
                        1,1,0,1,1,0,2,2,
                        2,1),
            save.context = c("CpG","CHG","CHH"),
            read.context="none",
            save.folder = "meth_read")

# Some SAM files were not fully processed it seems, so the process fails. Files are broken up to find which do not allow for the process to continue. Will come back to these after the deadline.
"0dph_rrbsDB.3" = processBismarkAln(file.list3,
            sample.id=list(), 
            assembly="omyk1.1A",
            treatment=c(),
            save.context = c("CpG","CHG","CHH"),
            read.context="none",
            save.folder = "meth_read")
```

```{r CpG methylation pt1, echo=FALSE, results=FALSE}
# Read in *CpG.txt files via methRead()
file.exists("C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/32_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/40_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/41_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/43_A_CpG.txt", 
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/46_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/47_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/50_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/54_A_CpG.txt", 
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/59_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/62_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/64_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/76_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/77_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/78_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/79_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/81_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/82_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/85_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/86_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/90_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/95_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/96_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/99_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/100_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/101_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/102_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/108_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/109_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/117_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/118_A_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/120_B_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/121_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/124_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/125_C_CpG.txt",
            "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/127_A_CpG.txt")

list_supp = list("C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/32_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/40_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/46_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/47_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/50_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/59_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/64_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/76_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/77_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/78_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/79_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/81_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/85_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/86_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/90_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/95_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/96_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/99_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/100_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/109_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/121_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/124_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/125_C_CpG.txt")
```

```{r CpG methylation pt2}
myobj_supp = methRead(list_supp,
         sample.id=list("32_Chol+","40_Chol+","46_Chol","47_Chol+","50_Chol","59_Chol+","64_Chol+","76_Chol",
                        "77_Chol","78_Chol+","79_Chol","81_Chol+","85_Chol","86_Chol+","90_Chol","95_Chol",
                        "96_Chol+","99_Chol+","100_Chol+","109_Chol+","121_Chol+","124_Chol+","125_Chol+"),
          assembly="omyk1.1A",
          treatment=c(1,1,0,1,0,1,1,0,
                      0,1,0,1,0,1,0,0,
                      1,1,1,1,1,1,1),
          context="CpG")
```

```{r CpG methylation pt3, echo=FALSE}
list_inadeq = list("C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/41_A_CpG.txt",
                   "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/43_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/46_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/50_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/54_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/62_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/76_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/77_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/79_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/82_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/85_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/90_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/95_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/101_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/102_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/108_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/109_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/117_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/118_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/127_A_CpG.txt")

myobj_inadeq = methRead(list_inadeq,
         sample.id=list("41_Chol-","43_Chol-","46_Chol","50_Chol","54_Chol-","62_Chol-","76_Chol","77_Chol",
                        "79_Chol","82_Chol-","85_Chol","90_Chol","95_Chol","101_Chol-","102_Chol-","108_Chol-",
                        "109_Chol","117_Chol-","118_Chol-","127_Chol-"),
          assembly="omyk1.1A",
          treatment=c(1,1,0,0,1,1,0,0,
                      0,1,0,0,0,1,1,1,
                      0,1,1,1),
          context="CpG")

list_supp_inadeq = list("C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/32_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/40_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/41_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/43_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/47_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/54_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/59_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/62_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/64_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/78_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/81_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/82_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/86_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/96_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/99_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/100_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/101_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/102_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/108_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/109_B_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/117_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/118_A_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/121_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/124_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/125_C_CpG.txt",
                "C:/Users/kfreij95/Documents/Objective 1 Choline Supplementation in Rainbow Trout in Vivo/meth_read/127_A_CpG.txt")

myobj_supp_inadeq = methRead(list_supp_inadeq,
         sample.id=list("32_Chol+","40_Chol+","41_Chol-","43_Chol-","47_Chol+","54_Chol-","59_Chol+","62_Chol-",
                        "64_Chol+","78_Chol+","81_Chol+","82_Chol-","86_Chol+","96_Chol+","99_Chol+","100_Chol+",
                        "101_Chol-","102_Chol-","108_Chol-","109_Chol+","117_Chol-","118_Chol-","121_Chol+","124_Chol+",
                        "125_Chol+","127_Chol-"),
          assembly="omyk1.1A",
          treatment=c(1,1,0,0,1,0,1,0,
                      1,1,1,0,1,1,1,1,
                      0,0,0,1,0,0,1,1,
                      1,0),
          context="CpG")
```
# Descriptive Statistics
```{r, echo=FALSE, eval=FALSE}
getMethylationStats(myobj_supp[[2]])
```
  
# Coverage and Differentially Methylated CpGs  
## Remove CpG sites with low coverage (below 10X coverage) and extremely high coverage > 99.9 percentile.  
```{r, message=FALSE}
#Supplemented vs Control
filtered.myobj.supp <- filterByCoverage(myobj_supp,lo.count=10,lo.perc=NULL,
                                   hi.count=NULL,hi.perc=99.9)
norm.myobj.supp <- normalizeCoverage(filtered.myobj.supp, method="median")
#Control vs Inadequate
filtered.myobj.inadeq <- filterByCoverage(myobj_inadeq,lo.count=10,lo.perc=NULL,
                                   hi.count=NULL,hi.perc=99.9)
norm.myobj.inadeq <- normalizeCoverage(filtered.myobj.inadeq, method="median")
#Supplemented vs Inadequate
filtered.myobj.supp.inadeq <- filterByCoverage(myobj_supp_inadeq,lo.count=10,lo.perc=NULL,
                                   hi.count=NULL,hi.perc=99.9)
norm.myobj.supp.inadeq <- normalizeCoverage(filtered.myobj.supp.inadeq, method="median")
```
  
# Comparative analysis  
## Merge samples
```{r, message=FALSE}
#Supplemented vs Control
meth.supp <- unite(norm.myobj.supp, destrand=FALSE, min.per.group = 2L)
#Control vs Inadequate
meth.inadeq <- unite(norm.myobj.inadeq, destrand=FALSE, min.per.group = 2L)
#Supplemented vs Inadequate
meth.supp.inadeq <- unite(norm.myobj.supp.inadeq, destrand=FALSE, min.per.group = 2L)

# Pool samples together to use Fisher exact test instead of beta regression
# pooled.meth <- pool(meth, sample.ids = c("C","A", "B"))
```
Number of bases per comparison:
- Supplemented-Control: 365955
- Control-Inadequate: 405470
- Supplemented-Inadequate: 555361

# PCA analysis
```{r, message=FALSE, echo=FALSE, results=FALSE}
#Source: https://yaaminiv.github.io/DML-Analysis-Part43/
plotColors <- rev(RColorBrewer::brewer.pal(5, "GnBu"))
plotCustomization <- data.frame("Group" = c("0"))

#Supplemented vs Control
allDataPCA.supp <- PCASamples(meth.supp, obj.return = TRUE)
summary(allDataPCA.supp)
#Look at summary statistics. The first PC explains 15.14% of variation, the second PC explains 7.75% of variation.
```
The PCA plots from methylKit weren't ideal, but the user above (Dr. Yaamini Venkataraman) has a great work around using the return object feature in the PCASamples to extract the data for our own customizable PCA plots.
```{r, message=FALSE, echo=FALSE}
fig.allDataPCA.supp <- ordiplot(allDataPCA.supp, choices=c(1,2), type = "none", display = "sites", 
                                xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object)
points(fig.allDataPCA.supp, col = c(rep(plotColors[2], times = 5), rep(plotColors[3], times = 5)), 
       pch = c(rep(16, times = 5), rep(17, times = 5)),"sites", cex = 2)
#Add each sample. Darker samples are control samples with lighter samples being samples from supplemented choline diet.
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 1, text = "PC 1 (15.14%)", line = 2)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 2, text = "PC 2 (7.75%)", line = 2)
legend("topright", legend = c("Control","Choline +"), pch = c(16,17), col = c(plotColors[2],plotColors[3]))
```
```{r, message=FALSE, echo=FALSE, results=FALSE}
#Control vs Inadequate
allDataPCA.inadeq <- PCASamples(meth.inadeq, obj.return = TRUE)
summary(allDataPCA.inadeq)
#Look at summary statistics. The first PC explains 18.37% of variation, the second PC explains 9.27% of variation.
```
```{r, message=FALSE, echo=FALSE}
fig.allDataPCA.inadeq <- ordiplot(allDataPCA.inadeq, choices=c(1,2), type = "none", display = "sites", 
                                xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object)
points(fig.allDataPCA.inadeq, col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), 
       pch = c(rep(16, times = 5), rep(18, times = 5)),"sites", cex = 2)
#Add each sample. Darker samples are choline deficient samples with lighter samples are control samples.
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")

axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 1, text = "PC 1 (18.37%)", line = 2)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 2, text = "PC 2 (9.27%)", line = 2)
legend("topright", legend = c("Control","Choline -"), pch = c(16,18), col = c(plotColors[2],plotColors[4]))
```
```{r, message=FALSE, echo=FALSE, results=FALSE}
#Supplemented vs Inadequate
allDataPCA.supp.inadeq <- PCASamples(meth.supp.inadeq, obj.return = TRUE)
summary(allDataPCA.supp.inadeq)
#Look at summary statistics. The first PC explains 12.08% of variation, the second PC explains 6.05% of variation.
```
```{r, message=FALSE, echo=FALSE}
fig.allDataPCA.supp.inadeq <- ordiplot(allDataPCA.supp.inadeq, choices=c(1,2), type = "none", display = "sites", 
                                xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object)
points(fig.allDataPCA.supp.inadeq, col = c(rep(plotColors[4], times = 5), rep(plotColors[3], times = 5)), 
       pch = c(rep(18, times = 5), rep(17, times = 5)),"sites", cex = 2)
#Add each sample. Darker samples are choline deficient samples with lighter samples are control samples.
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")

axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 1, text = "PC 1 (12.08%)", line = 2)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 2, text = "PC 2 (6.05%)", line = 2)
legend("topright", legend = c("Choline +","Choline -"), pch = c(18,17), col = c(plotColors[4],plotColors[3]))
```
  
## Get differentially methylated CpGs.
```{r, message=FALSE, warning=FALSE, results=FALSE}
#Supplemented vs Control
myDiff.supp <- calculateDiffMeth(meth.supp)
#Control vs Inadequate
myDiff.inadeq <- calculateDiffMeth(meth.inadeq)
#Supplemented vs Inadequate
myDiff.supp.inadeq <- calculateDiffMeth(meth.supp.inadeq)

# myDiff <- calculateDiffMeth(pooled.meth, num.cores = 4)
```
  
## Get hyper- and hypo-methylated CpGs.
### Get hypermethylated CpGs
```{r, message=FALSE, warning=FALSE, results=FALSE}
#Supplemented vs Control
myDiff25p.hyper.supp <- getMethylDiff(myDiff.supp,difference=20,qvalue=0.01,type="hyper")
myDiff25p.hyper.supp
#Control vs Inadequate
myDiff25p.hyper.inadeq <- getMethylDiff(myDiff.inadeq,difference=20,qvalue=0.01,type="hyper")
myDiff25p.hyper.inadeq
#Supplemented vs Inadequate
myDiff25p.hyper.supp.inadeq <- getMethylDiff(myDiff.supp.inadeq,difference=20,qvalue=0.05,type="hyper")
myDiff25p.hyper.supp.inadeq
```
  
Number of DMBs in Supplemented Comparison (q<0.05): 699
Number of DMBs in Supplemented Comparison (q<0.01): 286

Number of DMBs in Inadequate Comparison (q<0.05): 4993
Number of DMBs in Inadequate Comparison (q<0.01): 3433

Number of DMBs in Supplemented vs Inadequate (q<0.05): 
Number of DMBs in Supplemented vs Inadequate (q<0.01): 
  
### Get hypomethylated CpGs
```{r, message=FALSE, warning=FALSE, results=FALSE}
#Supplemented vs Control
myDiff25p.hypo.supp <- getMethylDiff(myDiff.supp,difference=20,qvalue=0.01,type="hypo")
myDiff25p.hypo.supp
#Control vs Inadequate
myDiff25p.hypo.inadeq <- getMethylDiff(myDiff.inadeq,difference=20,qvalue=0.01,type="hypo")
myDiff25p.hypo.inadeq
#Supplemented vs Inadequate
myDiff25p.hypo.supp.inadeq <- getMethylDiff(myDiff.supp.inadeq,difference=20,qvalue=0.05,type="hypo")
myDiff25p.hypo.supp.inadeq
```
  
Number of DMBs in Supplemented Comparison (q<0.05): 637
Number of DMBs in Supplemented Comparison (q<0.01): 266

Number of DMBs in Inadequate Comparison (q<0.05): 3830
Number of DMBs in Inadequate Comparison (q<0.01): 2721

Number of DMBs in Supplemented vs Inadequate (q<0.05): 
Number of DMBs in Supplemented vs Inadequate (q<0.01): 
  
### Get all Differentially Methylated Bases (DMBs)
```{r, message=FALSE, warning=FALSE, results=FALSE}
#Supplemented vs Control
myDiff25p.supp <- getMethylDiff(myDiff.supp,difference=20,qvalue=0.01)
myDiff25p.supp
#Control vs Inadequate
myDiff25p.inadeq <- getMethylDiff(myDiff.inadeq,difference=20,qvalue=0.01)
myDiff25p.inadeq
#Supplemented vs Inadequate
myDiff25p.supp.inadeq <- getMethylDiff(myDiff.supp.inadeq,difference=20,qvalue=0.05)
myDiff25p.supp.inadeq
```
  
Number of DMBs in Supplemented Comparison (q<0.05): 1336
Number of DMBs in Supplemented Comparison (q<0.01): 552

Number of DMBs in Inadequate Comparison (q<0.05): 8823
Number of DMBs in Inadequate Comparison (q<0.01): 6154

Number of DMBs in Supplemented vs Inadequate (q<0.05): 
Number of DMBs in Supplemented vs Inadequate (q<0.01): 
  
### Get data for all analyzed CpG sites
```{r, message=FALSE, warning=FALSE, results=FALSE}
#Supplemented vs Control 
myDiffall.supp <- getMethylDiff(myDiff.supp,difference=0,qvalue=1) 
myDiffall.supp
#Control vs Inadequate
myDiffall.inadeq <- getMethylDiff(myDiff.inadeq,difference=0,qvalue=1) 
myDiffall.inadeq
#Supplemented vs Inadequate 
myDiffall.supp.inadeq <- getMethylDiff(myDiff.supp.inadeq,difference=0,qvalue=1) 
myDiffall.supp.inadeq
```
