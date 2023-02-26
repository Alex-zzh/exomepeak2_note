## A function to extract exons grouped by unambiguous genes
exonsByiGenes <- function(txdb){
    exbg <- exonsBy(txdb, by = "gene") ##按照show所有exon，每个基因中的exon按照start从小到大顺序
    exbg <- exbg[elementNROWS(range(exbg)) == 1] ##只取每个基因所有exon最左端和最右端区间，保留
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
        return(unique_exons_gene)
    }
}
