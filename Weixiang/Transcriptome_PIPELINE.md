_PIPELINE by Weixiang Duan_ 
# Project: Transcriptome of <font color="red">_Littorina brevicula_</font> ![image][tmp]  
**RNASeq Analysis**   
**Last modified: 25th Oct, 2018**  

---  

## General Pipeline Step:  
  * [**Step_0: Backup file**](#0)  <font color="red">✓</font>
  * [**Step_1: FastQC**](#1)  <font color="red">✓</font>
  * [**Step_2: Trimmomatic**](#2)  <font color="red">✓</font>
  * [**Step_3: FastQC**](#3) <font color="red">✓</font>
  * [**Step_4: Trinity (2-5 days per population)**](#4) <font color="red">...</font>  
  * [**Step_5: CD-HIT-EST**](#5) 
  * [**Step_6: Quality assessment**](#6) 
    * (1) BUSCO 
    * (2) Transrate (2-3 days per population)
    * (3) ExN50  
  * [**Step_7: Salmon**](#7)
  * [**Step_8: TransDecoder**](#8)
    * (1) Extract the long open reading frames.
    * (2) Identify ORFs with homology to known proteins via blast or pfam searches (Search candidate peptides for homology).
        *  (a) Blastp Search
        *  (b) Pfam Search (2-3 days per population)
    * (3) Integrating the Blast and Pfam search results into coding region selection.
  * [**Step_9: ORTHOFINDER**](#9)
    * (1) OrthoFinder analysis (2-3 days)
    * (2) Prepare orthogroup-level count for each population
      * (a) Copy the gene-level counts file (quant.sf.genes) to the directory.
	    * (b) Prepare some files for the following analysis.
	    * (c) Use grep and awk to collect the counts of each OG in each sample, and then save them in one file (e.g. TL_counts.txt).
	    * (d) Rename the file (e.g. TL_counts.txt) for each sample.
	    * (e) Generate the final file, which shows the relationship between orthogroup and read counts of each sample.
    * (3) Merge OG2Counts files into one for the DESeq2 analysis in R.
  * [**Step_10: GO**](#10)
    * (1) INTERPRO (5-6 days per population)
    * (2) Prepare OG2GO file for each population
        * (a) Prepare EM_TL_pep2GO_ID2.txt and EM_TL_pep2OG_ID.txt
        * (b) Prepare EM_TL_pep2GO_part_1pep.txt
        * (c) Prepare TL_pep2GO_26617.txt
        * (d) Prepare EM_TL_OG2GO_26617.txt and EM_TL_pep2OG_part_1OG.txt
        * (e) Prepare TL_OG2GO_17835.txt
    * (3) Prepare OG2GO_Ontology file
  * [**Step_11: KEGG**](#11)  

---

## <span id="0">Step_0: Backup files. V  </span>
* Backup of RNA-seq data of *Littorina brevicula* (LB) were deposited in the following directory: (/mnt/sdd4t/DWX/LB_transcriptome_RawData_20180905), including both clean_data and raw_data.  

* Check the files completeness with MD5.txt files.
```MD5_script.sh``` script contains the following command lines:

```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome_RawData_20180905/raw_data/
dirlist=$(ls -d ${DIR}/*-LB-*)
cd $DIR
for dir in $dirlist
do
cd $dir
md5sum -c MD5*.txt >> ${DIR}/MD_log.txt
cd $DIR
done
```

* Build a new directory (/mnt/sdd4t/DWX/LB_transcriptome) for results deposits, and copy data into this directory.

```
cd /mnt/sdd4t/DWX
mkdir LB_transcriptome
cd /mnt/sdd4t/DWX/LB_transcriptome
mkdir 0_data
cd 0_data
cp /mnt/sdd4t/DWX/LB_transcriptome_RawData_20180905/raw_data/* /mnt/sdd4t/DWX/LB_transcriptome/0_data
```

* Rename file example:   
from "LYG\_LB\_25\_3\_**RRAS03252-V\_**1.fq.gz" to "LYG\_LB\_25\_3\_1.fq.gz"  
from "LYG\_LB\_25\_3\_**RRAS03252-V\_**2.fq.gz" to "LYG\_LB\_25\_3\_2.fq.gz"

## <span id="1">Step_1: FastQC</span> V  
* FastQC(Version 0.11.8) is used to assess the quality of NGS data.  
* Run the script ```Step1_FastQC.sh```    

```
cd /mnt/sdd4t/DWX/LB_transcriptome/1_FastQC
fromdos Step1_FastQC.sh
chmod +x Step1_FastQC.sh
./Step1_FastQC.sh        
```

* ```Step1_FastQC.sh``` script contains the following command lines:

```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/0_data
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/1_FastQC
cd /home/wj/LM_Analysis/FastQC
chmod 755 fastqc
cd $DIR
for sample in `ls *fq.gz`
do
/home/wj/LM_Analysis/FastQC/fastqc ${DIR}/$sample --outdir=$OUTDIR >> ${OUTDIR}/Step1_FastQC.log
done
```

## <span id="2">Step_2: Trimmomatic</span> V
* Trimmomatic is used to trim short reads and remove adapters.
* Run the script ```Step2_Trimmomatic.sh```

```
cd /mnt/sdd4t/DWX/LB_transcriptome/2_Trimmomatic
fromdos Step2_Trimmomatic.sh
chmod +x Step2_Trimmomatic.sh
./Step2_Trimmomatic.sh
```

* ```Step2_Trimmomatic.sh``` script contains the following command lines:

```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/0_data
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/2_Trimmomatic
cd $DIR
for sample in `ls *_1.fq.gz`
do
  base=$(basename $sample "_1.fq.gz")
  echo $base
  java -jar /home/wj/LM_Analysis/Trimmomatic/bin/trimmomatic.jar PE -threads 8 ${DIR}/${base}_1.fq.gz ${DIR}/${base}_2.fq.gz ${OUTDIR}/${base}_1.qc.fq.gz ${OUTDIR}/${base}_1.orphan.fq.gz ${OUTDIR}/${base}_2.qc.fq.gz ${OUTDIR}/${base}_2.orphan.fq.gz ILLUMINACLIP:/home/wj/LM_Analysis/Trimmomatic/adapters/Novogene_PE.fa:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25 >> ${OUTDIR}/Step2_trimmomatic.log
done
```

## <span id="3">Step_3: Run FastQC</span> V
* Run the script ```Step3_FastQC.sh```

```
cd /mnt/sdd4t/DWX/LB_transcriptome/3_FastQC
fromdos Step3_FastQC.sh
chmod +x Step3_FastQC.sh
./Step3_FastQC.sh        
```

* ```Step3_FastQC.sh``` script contains the following command lines:

```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/2_Trimmomatic
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/3_FastQC
cd /home/wj/LM_Analysis/FastQC
chmod 755 fastqc
cd $DIR
for sample in `ls *fq.gz`
do
/home/wj/LM_Analysis/FastQC/fastqc ${DIR}/$sample --outdir=$OUTDIR >> ${OUTDIR}/Step3_FastQC.log
done
```

## <span id="4">Step_4: Trinity</span>
### Trinity is used to perform de novo transcriptome assembly.
#### For XM population

```
cd /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity
fromdos Step4_Trinity_XM.sh
chmod +x Step4_Trinity_XM.sh
./Step4_Trinity_TL.sh        
```

* ```Step4_Trinity_TL.sh``` script contains the following command lines:

```
cd /mnt/sdd4t/DWX/LB_transcriptome/2_Trimmomatic
cat XM_LB_25_3_1.qc.fq.gz XM_LB_37_5_1.qc.fq.gz XM_LB_45_5_1.qc.fq.gz > /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/XM_LB_4inds_left.qc.fq.gz
cat XM_LB_25_3_2.qc.fq.gz XM_LB_37_5_2.qc.fq.gz XM_LB_45_5_2.qc.fq.gz > /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/XM_LB_4inds_right.qc.fq.gz
cd /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity
/home/wj/LM_Analysis/Trinity_v2.6.6/Trinity --seqType fq --left /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/XM_LB_4inds_left.qc.fq.gz --right /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_TL_4inds_right.qc.fq.gz --CPU 20 --max_memory 120G --bflyHeapSpaceMax 4G --SS_lib_type RF >> Step4_Trinity.TL.log
```

### For DT population

```
cd /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity
fromdos Step4_Trinity_DT.sh
chmod +x Step4_Trinity_DT.sh
./Step4_Trinity_XM.sh        
```

* ```/Step4_Trinity_XM.sh``` script contains the following command lines:

```
#!/bin/sh
cd /mnt/sdd4t/DWX/LB_transcriptome/2_Trimmomatic
cat XM-EM-25-2_1.qc.fq.gz XM-EM-37-10_1.qc.fq.gz XM-EM-45-9_1.qc.fq.gz XM-EM-52-5_1.qc.fq.gz > /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_XM_4inds_left.qc.fq.gz
cat XM-EM-25-2_2.qc.fq.gz XM-EM-37-10_2.qc.fq.gz XM-EM-45-9_2.qc.fq.gz XM-EM-52-5_2.qc.fq.gz > /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_XM_4inds_right.qc.fq.gz
cd /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity
/home/wj/LM_Analysis/Trinity_v2.6.6/Trinity --seqType fq --left /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_XM_4inds_left.qc.fq.gz --right /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_XM_4inds_right.qc.fq.gz --CPU 20 --max_memory 120G --bflyHeapSpaceMax 4G --SS_lib_type RF >> Step4_Trinity.XM.log
```

### For WZ population

```
cd /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity
fromdos Step4_Trinity_WZ.sh
chmod +x Step4_Trinity_WZ.sh
./Step4_Trinity_WZ.sh        
```

* this script contains the following command lines:

```
cd /mnt/sdd4t/DWX/LB_transcriptome/2_Trimmomatic
cat WZ-EM-25-3_1.qc.fq.gz WZ-EM-37-1_1.qc.fq.gz WZ-EM-45-1_1.qc.fq.gz WZ-EM-52-1_1.qc.fq.gz > /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_WZ_4inds_left.qc.fq.gz
cat WZ-EM-25-3_2.qc.fq.gz WZ-EM-37-1_2.qc.fq.gz WZ-EM-45-1_2.qc.fq.gz WZ-EM-52-1_2.qc.fq.gz > /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_WZ_4inds_right.qc.fq.gz
cd /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity
/home/wj/LM_Analysis/Trinity_v2.6.6/Trinity --seqType fq --left /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_WZ_4inds_left.qc.fq.gz --right /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_WZ_4inds_right.qc.fq.gz --CPU 20 --max_memory 120G --bflyHeapSpaceMax 4G --SS_lib_type RF >> Step4_Trinity.WZ.log
```

## <span id="5">Step_5: CD-HIT-EST</span>
### CD-HIT-EST is used to reduce the redundancy of assembly.

```
cd /mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST
cp /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/Trinity_TL/trinity_out_dir/Trinity.fasta ./
mv Trinity.fasta Trinity_TL.fasta
cp /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/Trinity_XM/trinity_out_dir/Trinity.fasta ./
mv Trinity.fasta Trinity_XM.fasta
cp /mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/Trinity_WZ//trinity_out_dir/Trinity.fasta ./
mv Trinity.fasta Trinity_WZ.fasta
```

### For each population (e.g. TL)

```
cd-hit-est -i Trinity_TL.fasta -o Trinity_TL_clustered_0.95.fa -c 0.95 -T 20 -M 120000 -n 10 -p 1 -g 1 -d 40 1>/mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.log 2>/mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.err
sed -i 's/TRINITY/TRINITY_TL/g' Trinity_TL_clustered_0.95.fa      # Change the name of transcripts
```

## <span id="6">Step_6: Quality assessment</span>

```
cd /mnt/sdd4t/DWX/LB_transcriptome/6_Quality
```

### (1) BUSCO
#### For each population (e.g. TL)

```
cd /mnt/sdd4t/DWX/LB_transcriptome/6_Quality/TL_BUSCO
run_BUSCO.py -i /mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.fa -o EM.TL_busco_metazoa -l /home/wj/LM_Analysis/BUSCO_v3/metazoa_odb9 -m transcriptome  --cpu 20
```

### (2) Transrate (2-3 days）  

#### For each population (e.g. TL)

```
cd /mnt/sdd4t/DWX/LB_transcriptome/6_Quality/TL_Transrate
sed 's_|_-_g' /mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.fa > /mnt/sdd4t/DWX/LB_transcriptome/6_Quality/TL_Transrate/EM.TL_fixed.fasta
transrate --assembly=/mnt/sdd4t/DWX/LB_transcriptome/6_Quality/TL_Transrate/EM.TL_fixed.fasta --left=/mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_TL_4inds_left.qc.fq.gz --right=/mnt/sdd4t/DWX/LB_transcriptome/4_Trinity/EM_TL_4inds_right.qc.fq.gz --output=/mnt/sdd4t/DWX/LB_transcriptome/6_Quality/TL_Transrate/6ty_v1 --threads 20
```

#### For each population (e.g. TL)
##### (a) Estimate transcript abundance with Salmon. Here, we can use the output files from Step_7 (/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL/quasi_index_k31_AllQuantSf/)to conduct this analysis.

```
cd /mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL/quasi_index_k31_AllQuantSf
cp TL-EM-25-5_quant.sf TL-EM-37-6_quant.sf TL-EM-45-11_quant.sf TL-EM-52-6_quant.sf /mnt/sdd4t/DWX/LB_transcriptome/6_Quality/TL_ExN50/
```

##### (b) Build transcript matrices.

```
cd /mnt/sdd4t/DWX/LB_transcriptome/6_Quality/TL_ExN50
/home/wj/LM_Analysis/Trinity_v2.6.6/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --out_prefix TL.salmon TL-EM-25-5_quant.sf TL-EM-37-6_quant.sf TL-EM-45-11_quant.sf TL-EM-52-6_quant.sf
```

##### (c) Count numbers of expressed transcripts.

```
/home/wj/LM_Analysis/Trinity_v2.6.6/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl TL.salmon.isoform.TPM.not_cross_norm | tee TL.salmon.isoform.TPM.not_cross_norm_by_min_TPM
less TL.salmon.isoform.TPM.not_cross_norm_by_min_TPM
```

##### (d) Contig Ex90N50 statistic and Ex90 transcript count.

```
/home/wj/LM_Analysis/Trinity_v2.6.6/util/misc/contig_ExN50_statistic.pl TL.salmon.isoform.TMM.EXPR.matrix /mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.fa | tee EM.TL_ExN50.stats
/home/wj/LM_Analysis/Trinity_v2.6.6/util/misc/plot_ExN50_statistic.Rscript  EM.TL_ExN50.stats
less EM.TL_ExN50.stats
```

## <span id="7">Step_7: Salmon</span>
### Use the following script, we can get transcript-level counts and gene-level-counts. This is different from our previous script which can only generate transcript-level counts.
#### For each population (e.g. TL)

```
cd /mnt/sdd4t/DWX/LB_transcriptome/7_Salmon
fromdos Step7_Salmon_TL_1_quasi_index_k31.sh
chmod +x Step7_Salmon_TL_1_quasi_index_k31.sh
./Step7_Salmon_TL_1_quasi_index_k31.sh        
```

* this script contains the following command lines:


```
#!/bin/sh
cd /mnt/sdd4t/DWX/LB_transcriptome/7_Salmon
# Build reference
/home/wj/LM_Analysis/Trinity_v2.6.6/util/align_and_estimate_abundance.pl --transcripts /mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.fa --est_method salmon --trinity_mode --prep_reference --salmon_idx_type quasi
# Perform analysis
DIR=/mnt/sdd4t/DWX/LB_transcriptome/2_Trimmomatic
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon
cd $OUTDIR
for sample in `ls $DIR/TL*_1.qc.fq.gz`
do
    name=$(basename $sample "_1.qc.fq.gz")
    echo $name
/home/wj/LM_Analysis/Trinity_v2.6.6/util/align_and_estimate_abundance.pl --transcripts /mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.fa \
--seqType fq --left ${DIR}/${name}_1.qc.fq.gz --right ${DIR}/${name}_2.qc.fq.gz \
--est_method salmon --SS_lib_type RF --trinity_mode --output_dir ${OUTDIR}/${name}_quant --thread_count 20
done
```

```
cd /mnt/sdd4t/DWX/LB_transcriptome/7_Salmon
fromdos Step7_Salmon_TL_2_quasi_index_k31_AllQuantLog.sh
chmod +x Step7_Salmon_TL_2_quasi_index_k31_AllQuantLog.sh
./Step7_Salmon_TL_2_quasi_index_k31_AllQuantLog.sh        
```

* ```/Step7_Salmon_TL_2_quasi_index_k31_AllQuantLog.sh``` script contains the following command lines:

```
#!/bin/sh
# This script is used to extract the salmon_quant.log file from each *_quant directory and rename them.
DIR=/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL/quasi_index_k31_AllQuantLog
cd $DIR
for sample in `ls -d *_quant`
do        
        base=$(basename $sample "_quant")
        echo $base        
        cd ${base}_quant
        cd logs
        mv salmon_quant.log ${base}_salmon_quant.log
        cp ${base}_salmon_quant.log $OUTDIR
        cd $DIR
done
cd /mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL/quasi_index_k31_AllQuantLog
grep "Mapping rate" *.log >  TL_salmon_quant_index_k31_MappingRate.txt
```

```
cd /mnt/sdd4t/DWX/LB_transcriptome/7_Salmon
fromdos Step7_Salmon_TL_3_quasi_index_k31_AllQuantSf.sh
chmod +x Step7_Salmon_TL_3_quasi_index_k31_AllQuantSf.sh
./Step7_Salmon_TL_3_quasi_index_k31_AllQuantSf.sh        
```

* this script contains the following command lines:

```
#!/bin/sh
# This script is used to rename salmon_quant.log file from each *_quant directory and copy it to quasi_index_k31_AllQuantLog directory.
DIR=/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL/quasi_index_k31_AllQuantSf
cd $DIR
for sample in `ls -d *_quant`
do        
        base=$(basename $sample "_quant")
        echo $base      
        cd ${base}_quant
        mv quant.sf ${base}_quant.sf
        cp ${base}_quant.sf $OUTDIR
        cd $DIR
done
```


## <span id="8">Step_8: TransDecoder</span> (https://github.com/TransDecoder/TransDecoder/wiki)
### For each population (e.g. TL)
#### (1) Extract the long open reading frames.

```
cd /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/1_LongestORF
cp /mnt/sdd4t/DWX/LB_transcriptome/5_CD-HIT-EST/Trinity_TL_clustered_0.95.fa ./
TransDecoder.LongOrfs -t Trinity_TL.v1_clustered_0.95.fa
cd Trinity_TL_clustered_0.95.fa.transdecoder_dir
cp longest_orfs.pep ../
cd /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/1_LongestORF
mv longest_orfs.pep EM.TL_longest_orfs.pep
cp EM.TL_longest_orfs.pep /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/
```

#### (2) Identify ORFs with homology to known proteins via blast or pfam searches (Search candidate peptides for homology).
##### (a) Blastp Search
* I use DIAMOND software to perform BLASTP.
* build reference. Since I have done this, others do not need to do this.

```
cd /mnt/sdd4t/WJ/Uniref90
diamond makedb --in uniref90.fasta -d uniref90
```

###### run diamond blastp

```
cd /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/2_Blastp
diamond blastp -p 20 -b 20  -f 6 --sensitive --max-target-seqs 1 --evalue 0.00001 --max-hsps 1 -d /mnt/sdd4t/WJ/Uniref90/uniref90.dmnd -q /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/EM.TL_longest_orfs.pep -o EM.TL_blastp.fa > EM.TL_blastp.log
```

##### (b) Pfam Search (2-3 days)
##### Search the peptides for protein domains using Pfam
```
cd /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/3_Pfam
hmmscan --cpu 2 --domtblout EM.TL_pfam.domtblout /mnt/sdd4t/WJ/Pfam/Pfam-A.hmm /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/EM.TL_longest_orfs.pep > EM.TL_pfam.log      # cpu最多使用2个
### (3) Integrating the Blast and Pfam search results into coding region selection.
#### The outputs generated above can be leveraged by TransDecoder to ensure that those peptides with blast hits or domain hits are retained in the set of reported likely coding regions.
cd /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/1_LongestORF
cp /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/2_Blastp/EM.TL_blastp.fa ./
cp /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/3_Pfam/EM.TL_pfam.domtblout ./
TransDecoder.Predict -t Trinity_TL_clustered_0.95.fa --retain_pfam_hits EM.TL_pfam.domtblout --retain_blastp_hits EM.TL_blastp.fa
cp Trinity_TL.v1_clustered_0.95.fa.transdecoder* /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/4_Predict/
#### The final coding region predictions will now include both those regions that have sequence characteristics consistent with coding regions in addition to those that have demonstrated blast homology or pfam domain content.
```

## <span id="9">Step_9: ORTHOFINDER</span>
### (1) OrthoFinder analysis
#### (a) Prepare input files

```
cd /mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TLvsXMvsWZ
cp /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/4_Predict/*pep ./
mv Trinity_TL_clustered_0.95.fa.transdecoder.pep TL_pep.fa
mv Trinity_XM_clustered_0.95.fa.transdecoder.pep XM_pep.fa
mv Trinity_WZ_clustered_0.95.fa.transdecoder.pep WZ_pep.fa
### (b) Run orthofinder (2~3 days)
cd /mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER
python /home/wj/LM_Analysis/OrthoFinder-master/orthofinder/orthofinder.py -f TLvsXMvsWZ -t 20 -a 10      # It seems that the workstation can only use 9 CPUs in blast even though I set 20 CPUs.
```

### (2) Prepare orthogroup-level count for each population
#### For each population (e.g. TL)
##### (a) Copy the gene-level counts file (quant.sf.genes) to the directory (/mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TL).
cd /mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TL
fromdos Step9_OrthoFinder_TL_1.sh
chmod +x Step9_OrthoFinder_TL_1.sh
./Step9_OrthoFinder_TL_1.sh              

* this script contains the following command lines:

```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TL
cd $DIR
for sample in `ls -d *_quant`
do      
        base=$(basename $sample "_quant")
        echo $base        
        cd ${base}_quant
        mkdir ${OUTDIR}/${base}_quant
        cp quant.sf.genes ${OUTDIR}/${base}_quant
        cd $DIR
done
```

#### (b) Prepare some files for the following analysis.
* Using R, prepare OG_unique.txt file (it contains the shared orthogroup ID by all three populations) according to the output file (Orthogroups.csv) of orthofinder analysis.

```
setwd("F:/E.malaccana_transcriptom/EM_Transcriptome_DataAnalysis/20_ORTHOGROUP/TLvsXMvsWZ/temp")
Orthogroups <- read.csv("Orthogroups.csv", head=TRUE, na.strings=c("","NA"))   # 81,242
Orthogroups_shared <- na.omit(Orthogroups)    # delete rows with empty value
write.csv(Orthogroups_shared, file="Orthogroups_shared.csv") # 43,163
OG_unique <- data.frame(Orthogroups_shared[,1])
colnames(OG_unique) <- c("OG_ID")
write.csv(OG_unique, file="OG_unique.csv")   # 43,163
```

#### On Windows, prepare TL_Gene2OG_uique.txt file according to the output file (Orthogroups.csv) (Orthogroups.csv --> TL_pep.OG.txt --> TL_pep.OG.xlsx --> TL_Gene2OG_uique.txt)(see TL_pep.OG.xlsx in the directory F:\E.malaccana_transcriptom\EM_Transcriptome_DataAnalysis\20_ORTHOGROUP )
##### Note: From TL_pep.OG.txt to TL_pep.OG.xlsx (pep2OG sheet), I use Excel to read TL_pep.OG.txt, then each pep_ID will be in a unique cell. The one-2-many relation can be converted to one-2-one relation (as shown in pep2OG sheet of TL_pep.OG.xlsx).  

#### (c) Use grep and awk to collect the counts of each OG in each sample, and then save them in one file (e.g. TL_counts.txt).
cd /mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TL
fromdos Step9_OrthoFinder_TL_2.sh
chmod +x Step9_OrthoFinder_TL_2.sh
./Step9_OrthoFinder_TL_2.sh              
* this script contains the following command lines:


```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TL
OGlist=$(cat ${DIR}/OG_unique.txt)
filelist=$(ls -d ${DIR}/*_quant)
cd $DIR
for sample in $filelist
do
cd $sample

for OGID in $OGlist
do
grep "$OGID" ${DIR}/TL_Gene2OG_unique.txt > TL_Gene2OG_${OGID}.txt
awk '{print $1}' TL_Gene2OG_${OGID}.txt > TL_TS_${OGID}.txt
grep -wf TL_TS_${OGID}.txt quant.sf.genes > TL_Gene2Counts_${OGID}.txt
awk '{print $5}' TL_Gene2Counts_${OGID}.txt > TL_Counts_${OGID}.txt
awk '{sum += $1};END {print sum}' TL_Counts_${OGID}.txt >> TL_counts.txt
done

cd $DIR
done
```

### (d) Rename the file (TL_counts.txt) for each sample.
```
cd /mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TL
fromdos Step9_OrthoFinder_TL_3.sh
chmod +x Step9_OrthoFinder_TL_3.sh
./Step9_OrthoFinder_TL_3.sh 
```             
* this script contains the following command lines:

```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/9_ORTHOFINDER/TL
cd $DIR
for sample in `ls -d ${DIR}/*_quant`
do        
        base=$(basename $sample "_quant")
        echo $base        
        cd ${base}_quant
        mv TL_counts.txt ${base}_counts.txt
        cp ${base}_counts.txt $DIR
        cd $DIR
done
```

### (e) Generate the final file, which shows the relationship between orthogroup and read counts of each sample.  
### Before conducting the following command, rename each xxx_counts.txt (e.g. change TL-EM-25-3_counts.txt to TL_EM_T25_S03_counts.txt)

```
paste OG_unique.txt *_counts.txt > TL_OG2Counts.txt
```

## (3) Merge OG2Counts files into one for the DESeq2 analysis in R.
### On Windows, generate EM_WZvsXMvsTL_OG2Counts.csv with Excel by merging the three files (TL_OG2Counts.txt, XM_OG2Counts.txt, WZ_OG2Counts.txt) and add name for each sample.


## <span id="10">Step_10: GO</span>
## (1) INTERPRO (https://github.com/ebi-pf-team/interproscan/wiki)
### For each population (e.g. TL)
#### Simplify the name of peptide sequence using the following command.
```
cd /mnt/sdd4t/DWX/LB_transcriptome/10_GO
awk '/^>/{print ">TL_ORF" ++i; next}{print}' < /mnt/sdd4t/DWX/LB_transcriptome/8_TransDecoder/EM.TL_longest_orfs.pep > /mnt/sdd4t/DWX/LB_transcriptome/10_GO/EM.TL_longest_orfs2.pep
```
#### Since the program (INTERPRO) will report error when * occurs in the sequence (see EM.TL_longest_orfs2.pep), we delete this character using the following command.
```
sed 's/\*//g' EM.TL_longest_orfs2.pep > EM.TL_longest_orfs3.pep   
```

#### run interproscan  (-6 days on LINEworkstation）

```
/home/wj/LM_Analysis/InterPro/interproscan-5.30-69.0/interproscan.sh -iprlookup -pa KEGG -goterms -i /mnt/sdd4t/DWX/LB_transcriptome/10_GO/EM.TL_longest_orfs3.pep  -b /mnt/sdd4t/DWX/LB_transcriptome/10_GO/EM.TL_Interpro_output -T /mnt/sdd4t/DWX/LB_transcriptome/10_GO/temp > TL_Interpro.log
```

## (2) Prepare OG2GO file for each population
### For each population (e.g. TL)  
#### (a) Prepare EM_TL_pep2GO_ID.txt with Excel by extracting two columns (pep_ID and GO_ID) from the output file (EM.TL_Interpro_output.tsv) from interproscan analysis. Then remove pep_ID without GO_ID and duplicates (i.e. the same pep_ID and GO_ID) in EM_TL_pep2GO_ID.txt with Excel.
####     Prepare EM_TL_pep2GO_ID2.txt from EM_TL_pep2GO_ID.txt with Excel by changing one-2-many relation to one-2-one relation, and remove duplicates.
####     Prepare EM_TL_pep2OG_ID.txt from Gene2OG sheet in TL_pep.OG.xlsx(see Step_9) by replacing "m." with "TL_ORF"

#### (b) Prepare EM_TL_pep2GO_part_1pep.txt (i.e. unique pep_ID with both GO_ID and OG_ID) with the following commands on R-Studio:
```
setwd("F:/E.malaccana_transcriptom/EM_Transcriptome_DataAnalysis/18_GO/InterproScan/temp")   # Use your own temporary directory
EM_TL_pep2GO_ID2 <- read.table("EM_TL_pep2GO_ID2.txt", head=FALSE)     # 96182 rows; 42980 unique pep_ID
EM_TL_pep2OG_ID <- read.table("EM_TL_pep2OG_ID.txt", head=FALSE)     # 80332 rows; 80332 unique pep_ID
EM_TL_pep2GO_part <- EM_TL_pep2GO_ID2[(EM_TL_pep2GO_ID2$V1 %in% EM_TL_pep2OG_ID$V1),]     # 60355 rows; 26617 unique pep_ID
EM_TL_pep2GO_part_1pep <- as.data.frame(unique(EM_TL_pep2GO_part$V1))  # get unique pep_ID (26617)
write.csv(EM_TL_pep2GO_part_1pep, file="EM_TL_pep2GO_part_1pep.csv")
```

#### (c) Prepare TL_pep2GO_26617.txt (i.e. the relation between pep_ID and GO_ID; one-2-many)
cd /mnt/sdd4t/DWX/LB_transcriptome/10_GO/temp
fromdos Step10_GO_TL_1.sh
chmod +x Step10_GO_TL_1.sh
./Step10_GO_TL_1.sh        

# this script contains the following command lines:
```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/10_GO/temp/TL
peplist=$(cat ${DIR}/EM_TL_pep2GO_part_1pep.txt)
cd $DIR
for pepID in $peplist
do
grep "\b${pepID}\b" EM_TL_pep2GO_ID2.txt > ${pepID}.txt
awk -vORS=, '{ print $2 }' ${pepID}.txt | sed 's/,$/\n/' >> TL_GO_26617.txt
done
rm TL_ORF*.txt
```
paste EM_TL_pep2GO_part_1pep.txt TL_GO_26617.txt > TL_pep2GO_26617.txt

#### (d) Prepare EM_TL_OG2GO_26617.txt and EM_TL_pep2OG_part_1OG.txt with the following commands on R-Studio:
```
setwd("F:/E.malaccana_transcriptom/EM_Transcriptome_DataAnalysis/18_GO/InterproScan/temp")
EM_TL_pep2OG_ID <- read.table("EM_TL_pep2OG_ID.txt", head=FALSE)
EM_TL_pep2GO_26617 <- read.table("TL_pep2GO_26617.txt", head=FALSE)
EM_TL_pep2OG_part <- EM_TL_pep2OG_ID[(EM_TL_pep2OG_ID$V1 %in% EM_TL_pep2GO_26617$V1),]
write.csv(EM_TL_pep2OG_part, file="EM_TL_pep2OG_part.csv")
EM_TL_pep2OG_part_1OG <- as.data.frame(unique(EM_TL_pep2OG_part$V2))  # get unique OG_ID
write.csv(EM_TL_pep2OG_part_1OG, file="EM_TL_pep2OG_part_1OG.txt")
EM_TL_pep2OG_part <- EM_TL_pep2OG_part[order(EM_TL_pep2OG_part$V1),]
EM_TL_OG2pep2GO_26617 <- cbind(as.character(EM_TL_pep2OG_part$V2), as.character(EM_TL_pep2OG_part$V1), as.character(EM_TL_pep2GO_26617$V2))
write.csv(EM_TL_OG2pep2GO_26617, file="EM_TL_OG2pep2GO_26617.csv")
```
#### Prepare EM_TL_OG2pep2GO_26617.txt from EM_TL_OG2pep2GO_26617.csv by modifying format.
#### Prepare EM_TL_OG2GO_26617.txt from EM_TL_OG2pep2GO_26617.txt by removing the pep_ID column.
#### Modifying the format of EM_TL_pep2OG_part_1OG.txt

#### (e) Prepare TL_OG2GO_17835.txt (Since there are duplicated OG_ID in EM_TL_OG2GO_26617.txt, we need make the OG_ID unique based on EM_TL_pep2OG_part_1OG.txt)
cd /mnt/sdd4t/DWX/LB_transcriptome/10_GO/temp
fromdos Step10_GO_TL_2.sh
chmod +x Step10_GO_TL_2.sh
./Step10_GO_TL_2.sh        # this script contains the following command lines:
```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/10_GO/temp/TL
OGlist=$(cat ${DIR}/EM_TL_pep2OG_part_1OG.txt)
cd $DIR
for OGID in $OGlist
do
grep "\b${OGID}\b" EM_TL_OG2GO_26617.txt > ${OGID}.txt
awk -vORS=, '{ print $2 }' ${OGID}.txt | sed 's/,$/\n/' >> test.txt    # extract the second column and make them into a row with comma separating them.
done
sed 's/\,/ /g' test.txt > test2.txt   # change the comma to blank
ruby -e 'STDIN.readlines.each { |l| l.split(" ").uniq.each { |e| print "#{e}," }; print "\n" }' < test2.txt > test3.txt   # remove duplicates of each row and separate them with comma
sed 's/\,$//g' test3.txt > EM_TL_OG_17835.txt   # remove comma at the end of each row
cp EM_TL_OG_17835.txt /mnt/sdd4t/DWX/LB_transcriptome/10_GO/temp
paste EM_TL_pep2OG_part_1OG.txt EM_TL_OG_17835.txt > TL_OG2GO_17835.txt
```

## (3) Prepare OG2GO_Ontology file (one-2-many)
## These two websites help us to retrieve GO_Ontology based on GO_ID. (http://www.geneontology.org/faq/how-do-i-get-term-names-my-list-go-ids) (https://yeastmine.yeastgenome.org/yeastmine/bag.do).
#### Prepare uniqueGO.txt (i.e. unique GO_ID) by merging files (TL_OG2GO_17835.txt, XM_OG2GO_17835.txt, WZ_OG2GO_17835.txt) and remove duplicates. Then upload this file to the website (https://yeastmine.yeastgenome.org/yeastmine/bag.do).
#### Prepare uniqueGO2Ontology.txt by modifying the format of download file. Prepare three files (GO2BP_ID.txt, GO2CC_ID.txt, GO2MF_ID.txt) with uniqueGO2Ontology.txt
#### For each Ontology, prepare a file (e.g. EM_GO2BP_ID.txt) with GO2BP_ID.txt and OG_unique.txt (see Step_9-(2)-(b))


## <span id="11">Step_11: KEGG</span>
## (a) KEGG ontology annotation
### For each population (e.g. TL)
cd /mnt/sdd4t/DWX/LB_transcriptome/11_KEGG
cp /mnt/sdd4t/DWX/LB_transcriptome/10_GO/EM.TL_longest_orfs2.pep ./
## upload EM.TL_longest_orfs2.pep to the website (https://www.kegg.jp/ghostkoala/) for KEGG annotation.
## Then we can download the output EM.TL.ko (pep_ID-2-KO_ID)

## (b) Prepare TL_KO_unique.txt and TL_KO2TS.txt from TL_TS2KO.xlsx (F:\E.malaccana_transcriptom\EM_Transcriptome_DataAnalysis\19_KEGG\TS2KO\TL)

## (c) Prepare TL_KO2Count.txt
### prepare TS2Count files
cd /mnt/sdd4t/DWX/LB_transcriptome/11_KEGG/TL
fromdos Step11_KEGG_TL_1.sh
chmod +x Step11_KEGG_TL_1.sh
./Step11_KEGG_TL_1.sh           

# this script contains the following command lines:
```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/7_Salmon/TL
OUTDIR=/mnt/sdd4t/DWX/LB_transcriptome/11_KEGG/TL/temp
cd $DIR
for sample in `ls -d *_quant`
do      
        base=$(basename $sample "_quant")
        echo $base        
        cd ${base}_quant
        mkdir ${OUTDIR}/${base}_quant
        cp quant.sf ${OUTDIR}/${base}_quant
        cd $DIR
done
```
### Change the directory name（TL-EM-25-3_quant ----> TL-EM-25-S03_quant）

### prepare TL_counts.txt for each sample
cd /mnt/sdd4t/DWX/LB_transcriptome/11_KEGG/TL
fromdos Step11_KEGG_TL_2.sh
chmod +x Step11_KEGG_TL_2.sh
./Step11_KEGG_TL_2.sh                 # this script contains the following command lines:           
```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/11_KEGG/TL/temp
KOlist=$(cat ${DIR}/TL_KO_unique.txt)
filelist=$(ls -d ${DIR}/*_quant)
cd $DIR
for sample in $filelist
do
cd $sample
for KOID in $KOlist
do
grep "$KOID" ${DIR}/TL_KO2TS.txt | awk '{print $2}' > TL_TS_${KOID}.txt
grep -wf TL_TS_${KOID}.txt quant.sf | awk '{print $5}' | awk '{sum += $1};END {print sum}' >> TL_counts.txt
done
rm TL_TS_K*.txt
cd $DIR
done
```

### prepare TL_KO2Count.txt by merging all individual sample
cd /mnt/sdd4t/DWX/LB_transcriptome/11_KEGG/TL
fromdos Step11_KEGG_TL_3.sh
chmod +x Step11_KEGG_TL_3.sh
./Step11_KEGG_TL_3.sh                 

# this script contains the following command lines:    
```
#!/bin/sh
DIR=/mnt/sdd4t/DWX/LB_transcriptome/11_KEGG/TL/temp

cd $DIR

for sample in `ls -d ${DIR}/*_quant`
do        
        base=$(basename $sample "_quant")
        echo $base        
        cd ${base}_quant
        mv TL_counts.txt ${base}_counts.txt
        cp ${base}_counts.txt $DIR
        cd $DIR
done

paste TL_KO_unique.txt *counts.txt > TL_KO2Count.txt
```

[tmp]:data:image/jpeg;base64,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
