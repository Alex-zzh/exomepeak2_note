# motif：RRACH等序列
# region：函数exonsByiGenes得到的结果
# sequence：物种BSgenome

## Motif sampler
sampleSequence <- function(motif, region, sequence, fixed = FALSE, replace = FALSE){
    #require(BSgenome)
    #require(GenomicFeatures)
    stopifnot(is(region, "GRangesList")|is(region, "GRanges"))
    if(is(region, "GRangesList")) region <- unlist(region)
    region <- reduce(region)
    
    ##获得所有基因对应的exon的序列
    region_dnass <- getSeq(x=sequence,
                            names=seqnames(region),
                            start=start(region),
                            end=end(region),
                            strand=strand(region),
                            as.character=FALSE)
    
    indx <- paste0("reg_", seq_along(region)) 
    regions_GRL <- split(region, indx)
    regions_GRL <- regions_GRL[indx] ##让每个exon区域对应一个名字
    rm(indx)
    vmp <- vmatchPattern(motif, region_dnass, fixed = fixed) ##模式匹配，可以看到motif长度也可以人选
    rm(region_dnass)
    vmp_gr <- GRanges(seqnames = rep(names(regions_GRL), elementNROWS(vmp)), ranges = unlist(vmp))
    rm(vmp)
    motif_on_regions <- mapFromTranscripts(vmp_gr,regions_GRL) ##获得motif对应的染色质位置
    rm(vmp_gr, regions_GRL)
    mcols(motif_on_regions) = NULL ##附加列值都不要
    seqlengths(motif_on_regions) <- seqlengths(region) ##seqlengths设置成exon区域的长度
    return(motif_on_regions)
}

##最终结果（以motif3为RRACH为例）
# GRanges object with 494799 ranges and 0 metadata columns:
#                    seqnames    ranges strand
#                       <Rle> <IRanges>  <Rle>
#     scaffold_1   scaffold_1 2132-2136      +
#     scaffold_1   scaffold_1 2251-2255      +
#     scaffold_1   scaffold_1 2289-2293      +
#     scaffold_1   scaffold_1 2390-2394      +
#     scaffold_1   scaffold_1 4820-4824      +
#            ...          ...       ...    ...
#   scaffold_993 scaffold_993 9280-9284      -
#   scaffold_993 scaffold_993 9262-9266      -
#   scaffold_993 scaffold_993 9231-9235      -
#   scaffold_993 scaffold_993 9164-9168      -
#   scaffold_993 scaffold_993 9143-9147      -
#   -------
#   seqinfo: 934 sequences from an unspecified genome; no seqlengths