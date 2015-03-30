## form URLs as
## http://sequenceontology.org/browser/current_release/term/SO:0001578
## 
## SO:0001632 downstream_gene_variant (intergenic)
## SO:0001631 upstream_gene_variant (intergenic)
## SO:0001619 nc_transcript_variant (in a non-coding RNA gene)
## SO:0001578 stop_lost (elongated transcript)
## SO:0001574 splice_acceptor_variant
## SO:0001575 splice_donor_variant

annotateSO <- function(variantsVR) {
  gSO <- sequence_variant.gSOXP
  nodeDataDefaults(gSO, "varIdx") <- integer(0)

  ## SO:0001819 synonymous_variant
  mask <- variantsVR$CONSEQUENCE %in% "synonymous"
  if (any(mask))
    nodeData(gSO, "SO:0001819", "varIdx") <- list("SO:0001819"=unique(variantsVR$VCFIDX[mask]))

  ## SO:0001627 intron_variant
  mask <- variantsVR$LOCATION %in% "intron"
  if (any(mask))
    nodeData(gSO, "SO:0001627", "varIdx") <- list("SO:0001627"=unique(variantsVR$VCFIDX[mask]))

  ## SO:0001587 stop_gained (premature stop codon)
  mask <- variantsVR$CONSEQUENCE %in% "nonsense"
  if (any(mask))
    nodeData(gSO, "SO:0001587", "varIdx") <- list("SO:0001587"=unique(variantsVR$VCFIDX[mask]))
 
  ## SO:0001589 frameshift_variant
  mask <- variantsVR$CONSEQUENCE %in% "frameshift"
  if (any(mask))
    nodeData(gSO, "SO:0001589", "varIdx") <- list("SO:0001589"=unique(variantsVR$VCFIDX[mask]))

  ## SO:0001583 missense_variant (amino acid change)
  mask <- variantsVR$CONSEQUENCE %in% "nonsynonymous"
  if (any(mask))
    nodeData(gSO, "SO:0001583", "varIdx") <- list("SO:0001583"=unique(variantsVR$VCFIDX[mask]))

  ## SO:0001624 3_prime_UTR_variant
  mask <- variantsVR$LOCATION %in% "threeUTR"
  if (any(mask))
    nodeData(gSO, "SO:0001624", "varIdx") <- list("SO:0001624"=unique(variantsVR$VCFIDX[mask]))
  
  ## SO:0001623 5_prime_UTR_variant
  mask <- variantsVR$LOCATION %in% "fiveUTR"
  if (any(mask))
    nodeData(gSO, "SO:0001623", "varIdx") <- list("SO:0001623"=unique(variantsVR$VCFIDX[mask]))

  ## SO:0001820 inframe_indel
  mask <- variantsVR$TYPE %in% c("Insertion", "Deletion") & variantsVR$LOCATION %in% "coding" &
          width(variantsVR) == 3
  if (any(mask))
    nodeData(gSO, "SO:0001820", "varIdx") <- list("SO:0001820"=unique(variantsVR$VCFIDX[mask]))

  gSO
}
