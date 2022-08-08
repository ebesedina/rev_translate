#!/usr/bin/env Rscript

# Name: reverseTranslateVariant.R
# Author: DMP
# Description: 
#  From a protein variant input (i.e. KRAS G12) reverese translates to a
#  genomic position.

# imports -----------------------------------------------------------------

library(magrittr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Biostrings)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# params ------------------------------------------------------------------

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

# data --------------------------------------------------------------------

input = c(
  "CEP97",
  "A9"
)

## testing -

## A1CF T360 -> chr10 52575851 25575853
## KRAS G12 -> chr12 25398283 25398285
## PCNP T98 -> chr3 101304293 101304295
## CEP97 A9 -> chr3 101443545 101443547
GeneAliases <- fread("/g/strcombio/fsupek_data/users/ebesedina/dCdF/data/GeneAliases_geneNames.txt", sep="\t", quote="") #new

# script ------------------------------------------------------------------
#not working if aminoacid is spanning 2 exons
rev_translate = function(gene_symbol, pchanges) {
  
  genes(txdb, columns = "GENEID") -> gns
  transcripts(txdb, columns = c("GENEID", "TXID")) -> txs
  cds(txdb, columns = c("GENEID", "TXID", "CDSID")) -> cds
  
  geneid = GeneAliases[`Approved symbol`==gene_symbol & Status =="Approved"]$`NCBI Gene ID`
  
  if (is.na(geneid)) { #for example C17orf113
    mapIds(
      x = org.Hs.eg.db, 
      keys = gene_symbol, 
      "ENTREZID", 
      "SYMBOL",
      multiVals = "first") -> geneid
  }
  geneid %<>% as.numeric() 
  
  ## filter the huge datasets
  gns %<>% .[as.numeric(.$GENEID) == geneid]
  ## why is there nas here?
  na_filt = !is.na(as.numeric(txs$GENEID))
  gid_filt =  as.numeric(txs$GENEID) == geneid
  txs %<>% .[na_filt & gid_filt] 
  
  cds$GENEID %>% lapply(function(x){
    any(x == geneid)
  }) %>% unlist() -> gid_filt
  cds %<>% .[gid_filt] 
  
  txids = txs$TXID
  if (length(txids) >0) {
    posible_seqs = rbindlist(lapply(txids, function(i) {
      # print(i)
      tx_name = i
      strd_g = strand(cds) %>% as.character() %>% unique()
      sel_cds = cds$TXID %>% lapply(function(x){any(x %in% i)}) %>% unlist()
      if (all(is.null(sel_cds))){ # my part
        ## this might happen if there is no coding seqence in transcript
        warning("No coding seq in transcript.")
        return(NULL)
      }
      if (all(!sel_cds)){
        ## this might happen if there is no coding seqence in transcript
        warning("No coding seq in transcript.")
        return(NULL)
      }
      cds_in_tx = cds[sel_cds]
      seq_cds = getSeq(genome,cds_in_tx)
      stopifnot(is.null(names(seq_cds)))
      if (strd_g == "-"){
        ## i am not sure, potential fail.
        order_idx = length(seq_cds):1
      } else {
        order_idx = 1:length(seq_cds)
      }
      fseq = DNAString()
      for (j in order_idx){
        fseq = c(fseq, seq_cds[[j]] )  
      }
      fps = translate(fseq)
      # better if we put a force here..
      # tst = as.character(fps[length(fps)]) == "*"
      if(as.character(fps[length(fps)]) != "*") {
        warning("wrong stop")
        return(NULL)}
      if(as.character(fps[1]) != "M") {
        warning("wrong start")
        return(NULL)}
      # is this always always true?
      # stopifnot(as.character(fps[1]) == "M")
      posible_seqs = lapply(pchanges, function(pchange) {
        protein_aa = pchange %>% stringr::str_extract("[:alpha:]")
        protein_pos = pchange %>% 
          stringr::str_extract("[:digit:]+") %>% 
          as.numeric()
        
        
        if (protein_pos > length(fps) ){
          exp_ch = "OUTOFPROT"
        } else {
          fps[protein_pos] %>% as.character() -> exp_ch
        }
        if (exp_ch != protein_aa){
          warning("wrong txs, checking next")
          return(NULL)
        } else {
          tx_pos_5p = protein_pos*3-2
          tx_pos_3p = protein_pos*3
          reorder_cds = cds_in_tx[order_idx]
          lengths(reorder_cds) %>% cumsum() -> csums
          tx_pos_5p_id = which(tx_pos_5p <= csums)[1]
          tx_pos_3p_id = which(tx_pos_3p <= csums)[1]
          tx_pos_5p_id
          ### here i remove all the extra positions that are contained in the
          ### other cds.
          if (tx_pos_5p_id == 1){
            base_pos_5p = tx_pos_5p
          } else {
            base_pos_5p = tx_pos_5p - csums[tx_pos_5p_id-1] 
          }
          
          if (tx_pos_3p_id == 1){
            base_pos_3p = tx_pos_3p
          } else {
            base_pos_3p = tx_pos_3p - csums[tx_pos_3p_id-1] 
          }
          
          if (strd_g == "-"){
            ## re-check these 1s!
            ge_pos_5p = end(reorder_cds[tx_pos_5p_id]) - base_pos_5p + 1
            ge_pos_3p = end(reorder_cds[tx_pos_3p_id]) - base_pos_3p + 1
          } else {
            ge_pos_5p = start(reorder_cds[tx_pos_5p_id]) + base_pos_5p - 1
            ge_pos_3p = start(reorder_cds[tx_pos_3p_id]) + base_pos_3p - 1
          }
          ## NOTE that the resulting positions do not necessarly be 3ntp long!
          ## for instance in splice site!
          ge_pos = GRanges(
            seqnames = unique(as.character(seqnames(reorder_cds))),
            ranges  = IRanges(
              start = min(c(ge_pos_5p,ge_pos_3p)),
              end = max(c(ge_pos_5p, ge_pos_3p))
            ),
            strand = as.character(strd_g)
          )
          
          ge_pos$pchange = pchange
          ge_pos$TXID = tx_name
          return(ge_pos)
        }
        
      })
      posible_seqs = posible_seqs[!sapply(posible_seqs,is.null)]
      posible_seqs = do.call("c", posible_seqs)
      posible_seqs = as.data.frame(posible_seqs) %>% data.table()
    }))
    if(nrow(posible_seqs)>0) {
      posible_seqs[, Gene:= gene_symbol]
      # posible_seqs[, TXID:=NULL]
      posible_seqs[, GENEID := geneid]
      posible_seqs = unique(posible_seqs)
      return(posible_seqs[])
    }
    
    
  }
}

rev_translate_toDNA <- function(gene_symbol, pchanges) {
  return(tryCatch(rev_translate(gene_symbol, pchanges), error=function(e) NULL))
}
