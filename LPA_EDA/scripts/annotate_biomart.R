#!/usr/bin/env Rscript
# Uso: Rscript scripts/annotate_biomart.R dados_processados/gene_ids.txt 103
# Release explícito: 103 é contemporâneo a 2021, NÃO confirmação do estudo.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop('Forneça arquivo com IDs e release Ensembl investigado')
if (!requireNamespace('biomaRt', quietly=TRUE)) stop('Instale BiocManager e biomaRt; consultar README')
ids_original <- readLines(args[1]); ids <- unique(sub('\\.[0-9]+$', '', ids_original))
release <- as.integer(args[2]); dir.create('anotacoes', showWarnings=FALSE)
keyfile <- tempfile(); writeLines(sort(ids), keyfile)
key <- unname(tools::md5sum(keyfile)); unlink(keyfile)
cache <- sprintf('anotacoes/cache_release%s_%s.rds',release,key)
if (file.exists(cache)) {
  result <- readRDS(cache)
} else {
  mart <- biomaRt::useEnsembl(biomart='genes',dataset='hsapiens_gene_ensembl',version=release)
  attrs <- c('ensembl_gene_id','external_gene_name','gene_biotype','description','chromosome_name','start_position','end_position','strand')
  pieces <- split(ids, ceiling(seq_along(ids)/500))
  result <- do.call(rbind,lapply(pieces,function(chunk) biomaRt::getBM(attributes=attrs,filters='ensembl_gene_id',values=chunk,mart=mart)))
  result$ensembl_release <- release
  datasets <- biomaRt::listDatasets(mart)
  result$genome_assembly <- datasets$version[match('hsapiens_gene_ensembl',datasets$dataset)]
  saveRDS(result,cache)
}
# Left join preserva não mapeados e relações um-para-muitos, sem escolher símbolos.
full <- merge(data.frame(original_ensembl_id=ids_original,ensembl_gene_id=sub('\\.[0-9]+$','',ids_original)),result,by='ensembl_gene_id',all.x=TRUE,sort=FALSE)
full$mapping_status <- ifelse(is.na(full$gene_biotype),'unmapped_or_missing_biotype','mapped')
write.csv(full,'anotacoes/biomart.csv',row.names=FALSE,na='')
writeLines(c(paste('Processed UTC:',format(Sys.time(),tz='UTC')),paste('Ensembl release:',release),'Compatibilidade com anotação original: ainda não confirmada; verificar métodos suplementares.',capture.output(sessionInfo())),'anotacoes/provenance.txt')
