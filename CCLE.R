
library('optparse')

option_list <- list(
    make_option(c("-g", "--gene_name"), type = "character", default=FALSE,
              help="given a gene you want to explore"),
    make_option(c("-t", "--ccle_tpm"), type = "character", default=FALSE,
              help="TPM or other expression index matrix of CCLE expression, rows coreesponds to genes, and columns to cellline name (default should be cellLine_Tissue)"),
    make_option(c("-c", "--ccle_count"), type = "character", default=FALSE,
              help="count matrix of CCLE expression for DESeq2, rows coreesponds to genes, and columns to cellline name (default should be cellLine_Tissue)"),
    make_option(c("-e", "--ccle_cell_ann"), type = "character", default=FALSE,
              help="optional, the cell annotation file that originally download from CCLE website"),
    make_option(c("-a", "--gene_ann"), type = "character", default=FALSE,
              help="optional, the gene annotation file where the gene_name and gene_type column should be included"),
    make_option(c("-p", "--prefix"), type = "character", default=FALSE,
              help="the prefix of output")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (all(do.call(c, opt) == F)){
  system('Rscript CCLE.R -h')
  quit()
  }

gene_name = opt$gene_name
prefix = opt$prefix
ccle_tpm_path = opt$ccle_tpm
ccle_count_path = opt$ccle_count
gene_ann_path = opt$gene_ann
ccle_cell_ann_path = opt$ccle_cell_ann

library('ggplot2')
library('DESeq2')
library('fgsea')
library(ggrepel)
library(ppcor)
# # parameters
# gene_ann_path = NA
# ccle_tpm_path = '/data1/XBH_Data/shiny-server-data/CCLE_exp/rnaseq/CCLE_RNAseq_rsem_genes_tpm_20180929.symbol.txt'
# ccle_count_path = '/data1/XBH_Data/shiny-server-data/CCLE_exp/rnaseq/CCLE_RNAseq_genes_counts_20180929.merge_gene.csv'
# ccle_cell_ann_path = '/data1/XBH_Data/shiny-server-data/CCLE_exp/rnaseq/Cell_lines_annotations_20181226.txt'
# gene_name = 'NCR3LG1'
# prefix = './B7/'

if (is.na(gene_ann_path)){
  static = './static/'
  gene_ann_path = paste0(static, 'gene_annotation_GENCODE.v22.csv')
}

# read file
read_in = function(ccle_tpm_path, ccle_count_path, ccle_cell_ann_path, gene_ann_path){
  gene_ann_cnt = read.csv(gene_ann_path) # CSV file for gene annotation
  protein_coding_ann = subset(gene_ann_cnt, gene_type == 'protein_coding')

  ## read CCLE tpm file, count file, and cell line annotation file
  # the file should be \t delimited and the first column coresponding to rownames
  tpm_mat = read.csv(ccle_tpm_path, sep = '\t', row.names = 1, check.names = F)

  count_matrix = read.csv(ccle_count_path, sep = '\t', row.names = 1, check.names = F)

  cell_ann = read.csv(ccle_cell_ann_path, sep = '\t', row.names = 1, check.names = F)

  return(list('protein_coding'=protein_coding_ann,
    'TPM' = tpm_mat,
    'Count' = count_matrix,
    'Cell' = cell_ann))
}
print("+++++ read in files")
inputRead = read_in(ccle_tpm_path, ccle_count_path, ccle_cell_ann_path, gene_ann_path)

protein_coding_ann = inputRead[['protein_coding']]
tpm_mat = inputRead[['TPM']]
count_matrix_all = inputRead[['Count']]
cell_ann = inputRead[['Cell']]

rm(inputRead)

print("+++++ get top and bottom expression")
## take the expression of given gene_ann_cnt
gene_exp = log2(as.vector(as.matrix(tpm_mat[gene_name,]))+1)
names(gene_exp) = colnames(tpm_mat)

## top and bottom 15% cuttoff
top_cutoff = quantile(gene_exp, 0.85)
bottom_cutoff = quantile(gene_exp, 0.15)


## look into the gene expression distribution
pdf(paste0(prefix, '_distribution_ccle_', gene_name, '.pdf'), width = 4, height = 4)
hist(gene_exp, 
    xlab = paste0('log2(TPM+1) of ', gene_name), 
    main = 'Histogram of CCLE Expression')
### add lines
abline(v = c(top_cutoff, bottom_cutoff),
        col = c('red', 'blue'), lty = 2)

legend('topright', legend=c('Top 15%', 'Bottom 15%'), pch = "-", col = c('red', 'blue'),
    cex = .7, bty = 'n')
dev.off()

## output top and bottom cell lines
gene_exp_out = data.frame("cellLine" = as.vector(names(gene_exp)), "expression" = gene_exp)

gene_exp_out = merge(gene_exp_out, cell_ann[,c("Name", "Pathology", "Site_Primary", "Histology")], by.x = "cellLine", by.y = 0)

gene_exp_out = gene_exp_out[order(-gene_exp_out$expression),]
write.table(gene_exp_out, file = paste0(prefix, "_expression_ccle.csv"), quote = F, row.names = F)
rm(gene_exp_out)


## ====== do correlation between all protein genes to the given gene ===========
print("+++++ calculating correlation of all genes")

#### partial correlation, adjusted by tissue type
tissue_types = lapply(strsplit(colnames(tpm_mat), '\\_'), function(x){paste0(x[2:length(x)], collapse = '_')})
tissue_types = do.call(c,tissue_types)
tissue_type_model_matrix = model.matrix(~0+factor(tissue_types))

partial_corr = function(x){
  pres = pcor(cbind(t(tpm_mat[c(x, gene_name),]), tissue_type_model_matrix))
  return(c('corr'=pres$estimate[1,2], 'pval'=pres$p.value[1,2]))
}
pcorr_res = sapply(sort(rownames(tpm_mat)[rownames(tpm_mat) != gene_name]), partial_corr)
pcorr_res = t(pcorr_res)
print("+++++ export correlation result")
## write out correlation
write.csv(pcorr_res, file = paste0(prefix, '_corr_', gene_name, '_othergenes.csv'), quote = F)


## do differential expression for gene top and bottom 15%
print('++++++ prepare data for DESeq2')
top = names(gene_exp[gene_exp>top_cutoff])
bottom = names(gene_exp[gene_exp<bottom_cutoff])
count_matrix = cbind(count_matrix_all[,top],
                     count_matrix_all[,bottom])

# using tissue type as batch, top and bottom as comparison
tissue_types = lapply(strsplit(colnames(count_matrix), '\\_'), function(x){paste0(x[2:length(x)], collapse = '_')})
tissue_types = do.call(c,tissue_types)
cond = DataFrame('batch' = tissue_types,
                   'cond' = c(rep(paste0(gene_name, '_Top'), length(top)), 
                              rep(paste0(gene_name, '_Bottom'), length(bottom))))
rownames(cond) = colnames(count_matrix)

print("+++++ running DESeq2")
## do DESeq2
dds = DESeqDataSetFromMatrix(count_matrix,
                             cond, ~ batch + cond)
dds = DESeq(dds)
dds <- readRDS(paste0(prefix, '_ccle_', gene_name, '_deseq_top_bottom.rds'))

print("+++++ export differential expression")
## get diff exp genes
compare_name <- paste0('cond_', gene_name, '_Top_vs_', gene_name, '_Bottom')
if (!compare_name %in% resultsNames(dds)){
    high_vs_low_ccle <- as.data.frame(results(dds, contrast = list(compare_name)))
}else{
    compare_name <- paste0('cond_', gene_name, '_Bottom_vs_', gene_name, '_Top')
    high_vs_low_ccle <- as.data.frame(results(dds, contrast = list(compare_name)))
    high_vs_low_ccle$log2FoldChange <- -high_vs_low_ccle$log2FoldChange
    high_vs_low_ccle$stat <- -high_vs_low_ccle$stat
}

write.table(high_vs_low_ccle, file = paste0(prefix, '_', gene_name, '_Top_vs_Bottom_ccle_deseq2_res_table.csv',
                                            quote = F)



## vocanoplot
vocano_plot <- function(mat){
  plot_df = subset(mat, (!is.na(padj)) & (!is.na(log2FoldChange)))
  plot_df$diff = (abs(plot_df$log2FoldChange) > log2(1.5)) & (plot_df$padj < 0.01)
  plot_df$diff = factor(plot_df$diff, levels = c('TRUE', 'FALSE'))
  labels = rbind(head(plot_df[order(plot_df$stat),], 20), 
                 tail(plot_df[order(plot_df$stat),], 30))
  
  labels$gene = rownames(labels)
  g = ggplot(data = plot_df,
         aes(x = log2FoldChange, y = -log10(padj), colour = diff))+
    geom_point()+#ylim(0, 30)+
    theme_bw()+
    theme(axis.text.x=element_text(size=10, colour = "black"),
          axis.text.y=element_text(size=10, colour = "black"),
          panel.border = element_blank(),axis.line = element_line(colour = "black"),
          text=element_text(size=14, colour = "black"),
          legend.text=element_text(size=10),
          plot.title = element_text(hjust=0.5,vjust = 0.5,
                                    margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))+
    scale_colour_manual(values=c('TRUE'='red', 'FALSE'='grey'))+
    geom_text_repel(data = labels, aes(x = log2FoldChange, y = -log10(padj), label = gene), colour = 'black')
  print(g)
}

print("+++++ draw volcano plot")
# volcano plot
pdf(paste0(prefix, '_vocanoplot_all.pdf'), width = 10, height = 6)
vocano_plot(high_vs_low_ccle[rownames(high_vs_low_ccle) %in% as.vector(protein_coding_ann$gene_name),])
dev.off()