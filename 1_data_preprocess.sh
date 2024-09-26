#!/bin/bash
# This is a custom script to process ATAC and ACC-seq data
# Zeran Jia, jzr18@mails.tsinghua.edu.cn         


wd=${workingdir}

# first build a sample list with all sample prefix
#   and fastq file need to be like fresh_1_R1.fastq.gz fresh_2_R2.fastq.gz

list_sample=''

# specify input & output directory
# all fastq file should be stored in one directory
# different sample will have different subdirectory under output_dir

fastq_dir=''
output_dir=''



#define path to reference genome mm10

#a sudo code define path to where you save the code to find reference genome expect $Index
wd=$(echo ${workingdir}/ref )


Index= #bowtie2 index, file is too big for email attachment
blacklist=$(echo ${workingdir}/ref/mm10.blacklist.withJDB.sorted.bed )
filter1='random'
filter2='chrUn'
TSS_BED=$(echo ${workingdir}/ref//mm10.vM23.pc.lv12.tss.bed.gz)
TSS_extend=$(echo ${workingdir}/ref/mm10.vM23.pc.lv12.tss.ext2k.bed)
CHROMSIZES=$(echo ${workingdir}/ref/mm10.chrom.sizes)



# major script

for sample in $(cat ${list_sample})
do

####  pre- process

align_dir=$(echo ${output_dir}/${sample}/align )
qc_dir=$(echo  ${output_dir}/${sample}/qc )
mkdir -p  ${align_dir}
mkdir -p  ${qc_dir}

##  detect  adaptor
# use python script in Kundaje pipeline, which was from https://github.com/nboley/GGR_code

python3 \
detect_adapter.py \
${fastq_dir}/${sample}_R1.fastq.gz \
> ${qc_dir}/adaptor1.log

python3 \
detect_adapter.py \
${fastq_dir}/${sample}_R2.fastq.gz \
> ${qc_dir}/adaptor2.log

adaptor_seq1=$(cat ${qc_dir}/adaptor1.log |sed -n 9p |cut -f 3 )
adaptor_seq2=$(cat ${qc_dir}/adaptor2.log |sed -n 9p |cut -f 3 )

## 
# cutadapt

cutadapt \
-j 3 -m 20 -e 0.1 -O 3 \
-q 20 --quality-base=33 \
-a ${adaptor_seq1} \
-A ${adaptor_seq2} \
-o ${align_dir}/${sample}_R1.trimmed.fastq.gz \
-p ${align_dir}/${sample}_R2.trimmed.fastq.gz \
${fastq_dir}/${sample}_R1.fastq.gz \
${fastq_dir}/${sample}_R2.fastq.gz \
> ${qc_dir}/${sample}_cutadapt_report.txt

## fastqc

fastqc \
-t 6 \
${fastq_dir}/${sample}_R1.fastq.gz \
${fastq_dir}/${sample}_R2.fastq.gz \
-o ${qc_dir} \
-d ${qc_dir}/


#### alignment
##   memory require : one  thread for bowtie2 is 5G

# about the alignment parameter -t -q -N1 -L 25
# could use --local or --end-to-end to replace

bowtie2 \
-X2000 \
--mm \
-t -q -N1 -L 25 --no-mixed --no-discordant \
--threads 8 \
-x ${Index} \
-1 ${align_dir}/${sample}_R1.trimmed.fastq.gz \
-2 ${align_dir}/${sample}_R2.trimmed.fastq.gz \
2>${qc_dir}/${sample}_bowtie2.log |\
samtools view -@ 8 -Su /dev/stdin |\
samtools sort -@ 8 -m 4G - > ${align_dir}/${sample}.bam

samtools index -@ 8 ${align_dir}/${sample}.bam

#### post-align

## filter 1/3    
# (use MAPQ instead of processing uniquely mapped reads;  uniquely mapping rarely mentioned today )
# flag: filter 1804=1024+512+256+8+4 ; get 2
# MAPQ > 30
# sort by name 
samtools view -F 1804 -f 2 -q 30 -@ 8 -u ${align_dir}/${sample}.bam |\
samtools sort -@ 8 -m 4G -n /dev/stdin -o ${align_dir}/${sample}.tmp.filt.bam

## filter 2/3
# fix mate info of name sorted bam
# sort by coordinate again 
samtools fixmate -@ 8 -r ${align_dir}/${sample}.tmp.filt.bam ${align_dir}/${sample}.tmp.fixmate.bam
samtools view -F 1804 -f 2 -@ 8 -u ${align_dir}/${sample}.tmp.fixmate.bam |\
samtools sort -@ 8 -m 4G /dev/stdin -o ${align_dir}/${sample}.filter.bam

# just need filter.bam for next ..
rm ${align_dir}/${sample}.tmp.filt.bam
rm ${align_dir}/${sample}.tmp.fixmate.bam

## filter 3/3
# picard mark duplicates (not remove) 
# use samtools view -F 1024(in 1804) to filter, better than picard ?

# had better specify a java temp path for Markduplicates 
#     or it might cause error when the default system path is full
# use the lastest version picard

java_temp=$(echo ${output_dir}/java_temp )
mkdir -p ${java_temp}

java -Xmx50G -Djava.io.tmpdir=${java_temp} \
-jar /conglilab/shared/applications/picard_2.20.4/picard.jar MarkDuplicates \
INPUT=${align_dir}/${sample}.filter.bam \
OUTPUT=${align_dir}/${sample}.dupmark.bam \
METRICS_FILE=${qc_dir}/${sample}.dup.qc \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=false

samtools view -F 1804 -f 2 -@ 8 -b -u ${align_dir}/${sample}.dupmark.bam |\
samtools sort -@ 8 -m 4G /dev/stdin -o ${align_dir}/${sample}.nodup.bam

samtools index -@ 8 ${align_dir}/${sample}.nodup.bam

## add one more step to filter blacklist
bedtools intersect -v -abam ${align_dir}/${sample}.nodup.bam -b ${blacklist} |\
samtools view -F 1804 -f 2 -@ 8 -S -h -b |\
samtools sort -@ 8 -m 4G  /dev/stdin -o  ${align_dir}/${sample}.nodup.clean.bam

samtools index -@ 8 ${align_dir}/${sample}.nodup.clean.bam

# flagstat ?
samtools flagstat -@ 8 ${align_dir}/${sample}.bam > ${qc_dir}/${sample}.flagstat
samtools flagstat -@ 8 ${align_dir}/${sample}.nodup.bam > ${qc_dir}/${sample}.nodup.flagstat

# library complexity
#    use dupmark.bam
echo "TotalPair,DictinctPair,OnePair,TwoPair,NRF=Distinct/Total,PBC1=OnePair/Distinct,PBC2=OnePair/TwoPair" >  ${qc_dir}/${sample}.pbc_qc.csv
bedtools bamtobed -i ${align_dir}/${sample}.dupmark.bam |\
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' |\
grep -v 'chrM' |sort |uniq -c |\
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} \
{m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; \
printf "%d,%d,%d,%d,%f,%f,%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' >>  ${qc_dir}/${sample}.pbc_qc.csv

# bigwig
# extend
bamCoverage \
-p 8 \
-e 250 \
-bs 100 \
--normalizeUsing RPKM \
-b ${align_dir}/${sample}.nodup.clean.bam \
-o ${align_dir}/${sample}.nodup.rpkm.bw

# remove temp files
#    or just hold for some time , may check for some possible usage
rm ${align_dir}/${sample}.filter.bam
rm ${align_dir}/${sample}.dupmark.bam

## nodup.bam to insert.bed
# filter with blacklist
# bam to bed ,and filter ,and shift 
# get insert.bed
# get insert.ext.bed , ext50 & ext250

bedtools intersect -v -abam ${align_dir}/${sample}.nodup.bam -b ${blacklist} |\
bedtools bamtobed -i - |gawk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' |\
grep -v chrM |grep -v ${filter1} |grep -v ${filter2} |\
gawk -F "\t" 'BEGIN{OFS=FS}{if($6=="+"){$2=$2+4} else if($6=="-"){$3=$3-5} print $0}' \
> ${align_dir}/${sample}.filter.shift.bed

cat ${align_dir}/${sample}.filter.shift.bed |gawk '{if($6=="+"){print $1"\t"$2-1"\t"$2;} else{print $1"\t"$3-1"\t"$3}}' \
> ${align_dir}/${sample}.insert.bed

cat ${align_dir}/${sample}.insert.bed |gawk '{print $1"\t"$3-50"\t"$3+50}' \
> ${align_dir}/${sample}.ext50.bed

cat ${align_dir}/${sample}.insert.bed |gawk '{print $1"\t"$3-250"\t"$3+250}' \
> ${align_dir}/${sample}.ext250.bed

#### add other qc like library complexity
#
## fragment_size histogram

java -Xmx30G -jar /conglilab/shared/applications/picard/picard.jar CollectInsertSizeMetrics \
I=${align_dir}/${sample}.nodup.clean.bam \
O=${qc_dir}/${sample}.nodup.fragsize.txt \
H=${qc_dir}/${sample}.nodup.fragsize.pdf \
VERBOSITY=ERROR QUIET=TRUE \
W=1000


## tss
# load conda env built for kundaje pipeline
#     dependency in this python script hasn't been independent yet
source activate bds_atac

OUTDIR=${qc_dir}/
OUTPREFIX=${sample}
FINALBAM=${align_dir}/${sample}.nodup.clean.bam

tss_plot='/conglilab/shared/projects/sc_omics/scripts/utils/tss_plot.py'

python ${tss_plot} \
--outdir $OUTDIR \
--outprefix $OUTPREFIX \
--tss $TSS_BED \
--finalbam $FINALBAM \
--chromsizes $CHROMSIZES

intersectBed -a ${TSS_extend} -b ${FINALBAM} |wc -l > ${qc_dir}/${sample}_reads_in_tss.txt
intersectBed -a ${TSS_extend} -b ${FINALBAM} -wa |sort -u |wc -l > ${qc_dir}/${sample}_reads_catched_tss.txt

source deactivate

## ATAC_qc

echo -e "Index,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag,Frag_mt,Frag_mt%,Frag_noMT,Alignment_rate,dup%,genomic_dup%,Tss_fold_enrichment,TSS_updown_2kb_reads,inTSS_ratio,TSS_updown_2kb_catched,catchTSS_ratio" \
> ${qc_dir}/${sample}.summary.csv

Reads=$(sambamba view ${align_dir}/${sample}.bam |wc -l )
Reads_m=$(sambamba view ${align_dir}/${sample}.bam |grep chrM |wc -l)
Reads_m_r=$(echo "scale=4;${Reads_m}/${Reads}" |bc)
Frag=$(sambamba view ${align_dir}/${sample}.nodup.bam |wc -l)
Frag_m=$(sambamba view ${align_dir}/${sample}.nodup.bam |grep chrM |wc -l)
Frag_m_r=$(echo "scale=4;${Frag_m}/${Frag}" |bc)
Frag_n=$(echo "${Frag}-${Frag_m}" |bc)
Align=$(cat ${qc_dir}/${sample}_bowtie2.log |grep "alignment rate" |gawk '{print $1}' )
Align_d=$(echo 0.01*${Align} |cut -d "%" -f 1 |bc)
Reads_n=$(printf "%.0f\n" `echo "${Reads}*${Align_d}-${Reads_m}" |bc`)
dup_r=$(echo "scale=4;1-${Frag}/(${Reads}*${Align_d})" |bc )
genomic_dup_r=$(echo "scale=4;(${Reads_n}-${Frag_n})/${Reads_n}" |bc )
TSS=$(printf $(cat ${qc_dir}/${sample}_tss-enrich.txt )"\n" )
inTSS=$(cat ${qc_dir}/${sample}_reads_in_tss.txt )
inTSS_r=$(echo "scale=4;${inTSS}/${Frag_n}" |bc)
catchTSS=$(cat ${qc_dir}/${sample}_reads_catched_tss.txt )
TSS_n=$(cat ${TSS_extend} |wc -l)
catchTSS_r=$(echo "scale=4;${catchTSS}/${TSS_n}" |bc)
echo -e ${sample}","${Reads}","${Reads_m}","${Reads_m_r}","${Reads_n}","${Frag}","${Frag_m}","${Frag_m_r}","${Frag_n}","${Align}","${dup_r}","${genomic_dup_r}","${TSS}","${inTSS}","${inTSS_r}","${catchTSS}","${catchTSS_r} \
>> ${qc_dir}/${sample}.summary.csv

done


## summary

cd  ${output_dir}
mkdir bam_file
mkdir bigwig_file
mkdir frag_size

echo -e "Index,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag,Frag_mt,Frag_mt%,Frag_noMT,Alignment_rate,dup%,Tss_fold_enrichment,TSS_updown2kb_reads,inTSS_ratio,TSS_updown_2kb_catched,catchTSS_ratio" \
> summary.csv

for sample in $(cat ${list_sample})
do

ln -s ${output_dir}/${sample}/align/${sample}.nodup.clean.bam* ./bam_file/
ln -s ${output_dir}/${sample}/align/${sample}.nodup.rpkm.bw ./bigwig_file/
ln -s ${output_dir}/${sample}/qc/${sample}.nodup.fragsize.pdf ./frag_size/

cat ${output_dir}/${sample}/qc/${sample}.summary.csv |tail -n 1 >> summary.csv
done

