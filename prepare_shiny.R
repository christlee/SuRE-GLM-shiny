library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(GenomicFeatures)
library(AnnotationDbi)
library(stringr)

gff_file = '/DATA/scratch/usr/c.leemans/data/tracks/hg19/gencode.v27lift37.annotation.gff3.gz'
rd_file = '/DATA/scratch/usr/c.leemans/projects/SuRE/SuRE_K562/promoter_triangle_input.RData'

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



ensembl_dt = max_dt[,list(ensembl_id = gsub('[.]*', '', gene_id),
                          transcript_id)]
setkey(ensembl_dt, 'ensembl_id')

gencode_dt = max_dt[, c('gene_id', 'transcript_id')]
setkey(gencode_dt, 'gene_id')

symbol_dt = max_dt[, c('gene_name', 'transcript_id')]
setkey(symbol_dt, 'gene_name')

transcript_dt = data.table(ensembl_tid = gsub('[.]*', '', names(tss_gr)),
                           transcript_id = names(tss_gr), key='ensembl_tid',
                           stringsAsFactors=F)


# lookup_list = list('symbol'=symbol_dt, 'ensembl_gene'=ensembl_dt,
#                    'ensembl_transcript'=transcript_dt, 'gencode'=gencode_dt)

lookup_list = list('symbol'=symbol_dt)
save(lookup_list, tss_gr, file=rd_file)
