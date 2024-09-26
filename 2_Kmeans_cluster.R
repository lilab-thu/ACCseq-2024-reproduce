#this sudo code is used to process ACC-seq data to identify hex-released peaks
#by Zeran Jia

#set environment ---------------------

library(dplyr)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(GenomicFeatures)

#set working directory
setwd('')
wd <- 'your working directory'


#merge replicate --------------------------

#this script is used to merge 2 replicate of each condition for further analysis

#merge bam file

for (condi in c('fresh','fix','hex')){

  #merge bam file
  cmd <- str_glue('samtools merge {condi}.nodup.clean.bam {condi}_1.nodup.clean.bam {condi}_2.nodup.cleanbam')
  system(cmd)
  
  cmd <- str_glue('samtools index -@ 8 {condi}.nodup.clean.bam')
  system(cmd)
}


#generate merged bigwig

for (condi in c('fresh','fix','hex')){
  
  cmd <- str_glue('bamCoverage -p 8 -bs 10 -e 250 -bl {blacklist} --normalizeUsing RPKM ',
  '-b {condi}.nodup.clean.bam -o {condi}.b10.nodup.rpkm.bw &')
  
}


#call peak for all samples-------------------------

for (condi in c('fresh','fix','hex')){
  
  #set path to merged bam
  merge_bam <- ''
  
  #define output directory
  subout_dir1 <- ''
  
  #merge bam file
  cmd <- str_glue('macs2 callpeak --verbose 3 --treatment {merge_bam} ',
                  '-g mm -B -q 0.05 --nomodel --keep-dup all -f BAM --outdir {subout_dir1} -n {condi}.nodup.noM')
  system(cmd)
}


#merge different peak together as input for following analysis

#define path to MACS2 narrowpeak
fresh_peak <- ''
fix_peak <- ''
hex_peak <- ''

cmd <- str_glue('multiIntersectBed -i {fresh_peak} {fix_peak} {hex_peak} > ROI_mti_combined.bed')



#following process:
#merge intersected peak to get a clean and clear input file
#bedtools sort -i ROI_mti_combined.bed >ROI_mti_combined.sorted.bed 
#bedtools merge -d 1000 -i ROI_mti_combined.sorted.bed >ROI_mti_combined.sorted.merge1k.bed



#use deeptools for data visualization and k-means clustering -------------------------

bw_list <-  c('path_to_fresh_bigwig',
              'path_to_fix_bigwig',
              'path_to_hex_bigwig')

ROI <- 'ROI_mti_combined.sorted.merge1k.bed'

oname <- str_glue('{wd}/res/ComputeMatrix.com_ROI_ab500.gz')
mname <- str_glue('{wd}/res/ComputeMatrix.com_ROI_ab500.txt')

cmd <- str_glue('computeMatrix scale-regions -S {bw_list[1]} {bw_list[2]} {bw_list[3]} -R {bed} -a 500 -b 500 ', 
                '-out {oname} -m 1500 --startLabel "peak start" --endLabel "peak end" ',
                '--outFileNameMatrix {mname} --samplesLabel fresh fix hex10')

system(cmd)


# out_dir = '/conglilab/shared/projects/ATAC_SDS/yinql_seq201117/analyze/210327_deeptools_cluster'

mname <- str_glue('{wd}/res/ComputeMatrix.com_ROI_ab500.gz')

bedname <- str_glue('{wd}/res/Heat500_k9.bed')
picname <- str_glue('{wd}/res/Heat500_k9.png')


cmd <- str_glue('plotHeatmap -m {mname} -o {picname} --outFileSortedRegions {bedname} --dpi 800 --alpha 0.8 --kmeans 9 ')

print(cmd)
system(cmd)

#!!!!!!important !!!!!!!!
#in this step, kmeans results may have some sligtly difference as difference of random seed
#as deeptools can not set consistent seed in different run

#in most cases, cluster 7 in c9 cluster is the most significant hex-released cluster, which will be C-I in cell paper figure
#cluster 1,2,3 after k9 kmeans are hex-consistant, which will be C-II, C-III and C-IV in cell paper figure


#split deeptools output bed file for following analysis

bedsum <- read_tsv(str_glue('{wd}/res/Heat500_k9.bed'))

colnames(bedsum) <- c('chr','start','end','name','score','strand','thickStart','thickEnd','itemRGB','blockCount','blockSize',
                      'blockStart','deepTools_group')

for (c in c(1:9)){
  clname <- str_glue('cluster_{c}')
  
  outname <- str_glue('{wd}/res/KmeansBed_km9_cl{c}.bed')
  
  bedsum %>% dplyr::filter(deepTools_group == clname) %>% dplyr::select(chr, start, end, deepTools_group) %>% 
    write_tsv(outname, col_names = F)
  
}


#GO analysis for hex-released peak
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(clusterProfiler)
  library(enrichplot)
  library(ChIPseeker)
  library(UpSetR)
  library(org.Mm.eg.db)
  library(ggplotify)
})

gtfanotation <- 'path to encode anotation gtf file'

gencodeV25<-makeTxDbFromGFF(file = "/conglilab/shared/genomics/mouse/gencode/encode/vm25/annotation/gencode.vM25.annotation.gtf")

released_peak <- str_glue('{wd}/res/KmeansBed_km9_cl7.bed')

released_peak_anno <- annotatePeak(released_peak, level="gene" , tssRegion=c(-3000, 3000), TxDb=gencodeV25,annoDb = "org.Mm.eg.db") %>% 
  as.data.frame()


# define cl7 genes

cl7_genes <- released_peak_anno %>% na.omit() %>% dplyr::pull(SYMBOL)
cl7_genes %>% length()

# run go
cl7_go_result <- enrichGO(
  gene = cl7_genes,
  keyType = "SYMBOL",
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  pvalueCutoff = 5e-2,
  pAdjustMethod = "BH",
  qvalueCutoff = 1e-1,
  minGSSize = 10,
  maxGSSize = 500,
)

cl7_go_result %>% as.data.frame() %>% write_csv(str_glue('{wd}/res/eGO.k9.c7.csv'), col_names = T)

# read predefined development realted pathways
dev_related_pathways <- read_csv("path/to/defined_terms.txt", col_names = FALSE) %>%
  pull(X1) %>%
  factor(levels = .)
dev_related_pathways


filter_go_result <- function(go_result) {
  go_result %>%
    mutate(Description = str_to_lower(Description)) %>%
    filter(Description %in% levels(dev_related_pathways)) %>%
    mutate(Description = factor(Description, levels = levels(dev_related_pathways))) %>%
    arrange(desc(Description))
}

wrap_char_count <- 30
cl7_go_result %<>%
  as.data.frame() %>%
  as.tibble() %>%
  filter_go_result()

transform_go_results <- function(go_result, name) {
  orig <- go_result
  tibble(
    name = name,
    id = orig$ID,
    term = str_wrap(orig$Description, width = wrap_char_count),
    padj_score = -log10(orig$pvalue),
    count = orig$Count
  )
}

plot_tbl <- cl7_go_result %>%
  transform_go_results(., "cl7")

plot_tbl %>%
  ggplot(.x, aes(x = padj_score, y = term)) +
  geom_col(width = 0.6, color = "black", aes(fill = count)) +
  scale_fill_gradient2(
    name = str_wrap("count", width = 15),
    low = "darkblue", high = "red", mid = "white",
    midpoint = 100,
    limits = c(0, 300),
    breaks = seq(0, 300, 100),
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, max_count)
  ) +
  ggtitle(.x$name %>% unique()) +
  theme_classic() +
  theme(
    aspect.ratio = 2,
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  )

















