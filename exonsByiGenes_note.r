## A function to extract exons grouped by unambiguous genes
exonsByiGenes <- function(txdb){
    exbg <- exonsBy(txdb, by = "gene") ##按照show所有exon，每个基因中的exon按照start从小到大顺序
    exbg <- exbg[elementNROWS(range(exbg)) == 1]
    fol <- findOverlaps(exbg) ##寻找基因区间之间的交集
    fol <- fol[queryHits(fol) != subjectHits(fol)] ##找出基因区间之间有交叉的基因
    ol_indx_M <- as.matrix(fol) ##取出所有有交集的基因区间对
    if (nrow(ol_indx_M) == 0) {
        return(reduce(exbg)) ##返回合并交集的基因exon区间
    }
    ##else部分整体要做的就是把有交叉的基因名都保存下来，所谓基因分区的名字，其他部分保持不变
    else {
        rm(fol)
        new_gene_names_temp <- names(exbg)
        new_gene_names_list <- split(new_gene_names_temp, seq_along(new_gene_names_temp)) 
        #R语言中的 seq_along() 函数用于生成一个与所传参数长度相同的序列。
        #这一步做的相当于把每个基因分成一类
        
        for (i in seq_len(nrow(ol_indx_M))) { #R 语言中的 seq_len() 函数用于生成从 1 到指定数字的序列。
        temp_i <- ol_indx_M[i, 1]
        new_gene_names_list[[temp_i]] <- c(new_gene_names_list[[temp_i]],
                                            new_gene_names_temp[ol_indx_M[i, 2]])
        } #这个循环相当于将有交叉的基因的对应交叉的基因两个基因名一同存下来，没有交叉的基因只存基因本身基因名
        
        rm(ol_indx_M, temp_i, new_gene_names_temp)
        new_gene_names_list <- lapply(new_gene_names_list, sort)
        new_gene_names <- vapply(new_gene_names_list, function(x) paste(x,
                                                                        collapse = ","), character(1))
        names(exbg) <- new_gene_names
        rm(new_gene_names, new_gene_names_list)
        rd_exons <- reduce(unlist(exbg), min.gapwidth = 0L) 
        #Ranges separated by a gap of at least min.gapwidth positions are not merged.
        
        fol <- findOverlaps(rd_exons, exbg) #相当于找到exon和每个基因的交叉关系
        split_indx <- rep(NA, length(rd_exons))
        split_indx[queryHits(fol)] <- names(exbg)[subjectHits(fol)] ##相当于给每个exon分配基因
        unique_exons_gene <- s plit(rd_exons, split_indx) #按照基因分配exon
        return(unique_exons_gene) ##最终得到的还是一个基因对应多个exon区间
    }
}

##最终结果
# GRangesList object of length 11395:
# $EWB00_000001
# GRanges object with 6 ranges and 0 metadata columns:
#         seqnames    ranges strand
#            <Rle> <IRanges>  <Rle>
#   [1] scaffold_1 2132-2441      +
#   [2] scaffold_1 2485-2539      +
#   [3] scaffold_1 4792-4978      +
#   [4] scaffold_1 6046-6156      +
#   [5] scaffold_1 6184-6288      +
#   [6] scaffold_1 6502-6641      +
#   -------
#   seqinfo: 934 sequences from an unspecified genome; no seqlengths

# $EWB00_000002
# GRanges object with 4 ranges and 0 metadata columns:
#         seqnames      ranges strand
#            <Rle>   <IRanges>  <Rle>
#   [1] scaffold_1   7206-7320      +
#   [2] scaffold_1   8939-9070      +
#   [3] scaffold_1 10259-10450      +
#   [4] scaffold_1 11608-11885      +
#   -------
#   seqinfo: 934 sequences from an unspecified genome; no seqlengths

# $EWB00_000003
# GRanges object with 7 ranges and 0 metadata columns:
#         seqnames      ranges strand
#            <Rle>   <IRanges>  <Rle>
#   [1] scaffold_1 19184-19700      -
#   [2] scaffold_1 20355-21213      -
#   [3] scaffold_1 22013-22108      -
#   [4] scaffold_1 27481-27581      -
#   [5] scaffold_1 27619-27771      -
#   [6] scaffold_1 27809-28108      -
#   [7] scaffold_1 30087-30296      -
#   -------
#   seqinfo: 934 sequences from an unspecified genome; no seqlengths

# ...
# <11392 more elements>