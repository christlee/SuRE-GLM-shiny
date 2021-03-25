library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(GenomicFeatures)
library(AnnotationDbi)
library(stringr)

gff_file = '/DATA/scratch/usr/c.leemans/data/tracks/hg19/gencode.v27lift37.annotation.gff3.gz'
tss_file = '/DATA/scratch/usr/c.leemans/projects/SuRE/SuRE_K562/tss_matrix.txt'
lookup_matrix = '/DATA/scratch/usr/c.leemans/projects/SuRE/SuRE_K562/lookup_matrix.txt'

fantom_file = paste0('/DATA/scratch/usr/c.leemans/projects/SuRE/tss_selection_gene_name.txt')
chrom_vec = paste0('chr', c(1:22, 'X', 'Y'))

txdb = makeTxDbFromGFF(gff_file)

tss_gr = promoters(txdb, columns=c('gene_id', 'tx_name'),
                   upstream=0, downstream=1,
                   filter=list(tx_chrom=chrom_vec))

tss_dt = fread(fantom_file, sep='\t',
               col.names=c('seqnames', 'pos', 'transcript_id',
                           'strand', 'gene_id', 'mean_expr', 'tissues_expr',
                           'gene_name'))

max_dt = tss_dt[, .SD[order(tissues_expr, mean_expr, decreasing=T)[1],],
                by='gene_id']

max_dt[, ensembl_id := gsub('[.]*', '', gene_id)]



fwrite(max_dt, file=lookup_matrix)
fwrite(as.data.table(tss_gr), file=tss_file)
