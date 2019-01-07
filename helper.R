# after ensembl query
# remove duplicated entries if any
# replace NA with original ensembl_gene_id
convert_mouse_ensembl_to_mgi_symbol <- function(ensembl_gene_id){
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  val <- ensembl_gene_id 
  out <- getBM(attributes=c('ensembl_gene_id', "mgi_symbol"),
    filters = 'ensembl_gene_id',
    values = val,
    mart = ensembl)

  unmatched <- genelist[! genelist %in% res[,1] ]
  unmatched <- cbind(unmatched, toupper(unmatched))
  colnames(unmatched) <- colnames(res)
  res <- rbind(res, unmatched)

  # resort the result to match original order
  idx <- match(genelist, res[,1])
  return(res[idx, ])

  return(out)
}

convert_mouse_gene_to_human_gene_symbol <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  #mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  res <- getLDS(attributes = c("mgi_symbol"),
          filters = "mgi_symbol",
          values = genelist,
          mart = mouse,
          attributesL = c("hgnc_symbol"),
          martL = human
          )
  unmatched <- genelist[! genelist %in% res[,1] ] %>% unlist() %>% as.character()
  unmatched <- cbind(unmatched, toupper(unmatched))
  colnames(unmatched) <- colnames(res)
  res <- rbind(res, unmatched)

  # resort the result to match original order
  idx <- match(genelist, res[,1])
  return(res[idx, ])
}

#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

create.brewer.color <- function(data, num=8, name="Set1")
{
	if(name=="naikai"){
		my_pallete <- c(
                        rgb(236,67,35, maxColorValue=255),
                        rgb(56,146,208, maxColorValue=255),
                        rgb(230,235,88, maxColorValue=255),
                        rgb(116,187,88, maxColorValue=255),
                        rgb(196,95,46, maxColorValue=255),
                        rgb(203,77,202, maxColorValue=255),
                        rgb(118,220,65, maxColorValue=255),
                        rgb(115,113,206, maxColorValue=255))
		return(create.manual.color(data, my_pallete))
	}else if(name=="naikai2"){
		my_pallete <- c(
                        rgb(10,10,10, maxColorValue=255),
                        rgb(230,235,88, maxColorValue=255),
                        rgb(118,220,65, maxColorValue=255),
                        rgb(203,77,202, maxColorValue=255),
                        rgb(196,95,46, maxColorValue=255),
                        rgb(116,187,88, maxColorValue=255),
                        rgb(236,67,35, maxColorValue=255),
                        rgb(56,146,208, maxColorValue=255),
                        rgb(115,113,206, maxColorValue=255))
		return(create.manual.color(data, my_pallete))
	}else if(name=="cellline"){
		my_pallete <- c(
                        "#89EB6B",
                        "#349B1B",
                        "#FF7265",
                        "#DA364D",
                        "#6099FF",
                        "#2256DB",
                        "#BEBADA",
                        "#00006D",
                        rgb(230,235,88, maxColorValue=255),
                        rgb(118,220,65, maxColorValue=255),
                        rgb(203,77,202, maxColorValue=255),
                        rgb(196,95,46, maxColorValue=255),
                        rgb(116,187,88, maxColorValue=255),
                        rgb(236,67,35, maxColorValue=255),
                        rgb(56,146,208, maxColorValue=255),
                        rgb(115,113,206, maxColorValue=255))
		return(create.manual.color(data, my_pallete))
	}else{
		groupCodes <- as.factor(data)
		colorCodes <- colorRampPalette(brewer.pal(num, name))(num)
		color.idx <- match(groupCodes, levels(groupCodes))
		label.color <- colorCodes[color.idx]
		return(label.color)
	}
}


#' Manually specify coloring for provided groups 
#'
#' Assign colors to data, if more unique data than provided colors, will impute the missing ones and fill them in
#' @param data data
#' @param group.color colors for each factorial group
#' @keywords color manual
#' @export
#' @examples
#' create.manual.color(c(1,2,3,1,2,3,2,3), c("red", "blue", "green"))

create.manual.color <- function(data, group.color)
{
	num_uniq_data <- length(unique(data))
	if(num_uniq_data > length(group.color)){
		require(RColorBrewer)
		group.color <- colorRampPalette(group.color)(num_uniq_data)
	}
	data <- as.numeric(factor(data))
	return(group.color[data])
}



#' Check rownames from data and auto detect which platform it is from
#'
#' Basically we want to convert those affy metrix array data and just
#' leave normal gene_symbols intact
#' @param data Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' detect_genecnt_platform("TP53")

detect_genecnt_platform <- function(data, method="mean"){
	require("annotate")
	require(magrittr)

	PROBES <- as.character(rownames(data))
	num.genes <- nrow(data)

	# for now just use simple nrow as check point
	if(num.genes == 22283){
  	  require(hgu133a.db) # 22283 probes
	  GPL="GPL96"
	  print ("Guess it's hgu133a")
	}else if(num.genes == 22645){
  	  require(hgu133b.db) # 22645 probes
	  GPL="GPL97"
	  print ("Guess it's hgu133b")
	}else if(num.genes == 54675 || num.genes == 44137){
  	  require(hgu133plus2.db) # 54675 probes
  	  # hard fix here for now, need to add mapped rownames to all these probe.annotation and then decide which one 
	  GPL="GPL570"
	  print ("Guess it's hgu133plus2")
	}else if(num.genes == 22277){
  	  require(hgu133a2.db) 
	  GPL="GPL571"
	  print ("Guess it's hgu133a2")
	}else{
	  print ("Guess its normal RNA-Seq, do nothing")
	  return(data)
	}

	if(method == "mean"){
	  select.fun <- function(x) rowMeans(x)
	}else if(method == "sd"){
	  select.fun <- function(x) rowSds(x)
	}else if(method == "max"){
	  select.fun <- function(x) rowMax(x)
	}else if(method == "min"){
	  select.fun <- function(x) rowMin(x)
	}else if(method == "iqr"){
	  select.fun <- function(x) apply(x, 1, IQR)
	}
	# There are two levels of mutliple mapping
	# one is that each probe_ids may have multiple mapped gene_symbols
	# For now, we just select the first gene_symbols if there are multiples
	gene.symbols <- affy_probe_to_gene_symbol(PROBES, GPL) %>%
	                filter(!duplicated(PROBEID)) %>%
	                dplyr::select(SYMBOL) %>%
	                unlist

	# Another place is that there may be multiple probes designed for each gene
	# Now we select the probes with the max(mean expression across samples)
	# Mayb fix it later by adding more different filter selections? Ray. 2015-10-30
	data.gene <- data %>% as.data.table %>%
	                      '['(, RowExp := select.fun(.SD)) %>%
	                      '['(, SYMBOL := gene.symbols) %>% group_by(., SYMBOL) %>%
	                      filter(., which.max(RowExp))

	# Remove rows without gene_symbols (NAs), then convert back to data.frame
	final.data <- data.gene %>%
	                  filter(!is.na(SYMBOL)) %>%
	                  setDF %>%
	                  set_rownames(.$SYMBOL) %>%
	                  '['(, !(colnames(.) %in% c("RowExp", "SYMBOL")))

	return(final.data)
}

#' A Cat Function
#'
#' This function allows you to read data with fread (faster)
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
# plot tSNE result 
myfread.table <- function(filepath, check.platform=T, header=T, sep="\t", detect.file.ext=T){
   require(data.table)
   require(tools)
   require(magrittr)

   ext <- file_ext(filepath)

   if(detect.file.ext){
      if (ext=="csv"){
        sep=","
      }else if (ext=="out" || ext=="tsv" || ext=="txt"){
          sep="\t"
      }else{
          warning("File format doesn't support, please try again")
          return(NULL)
      }
   }

   header <- read.table(filepath, nrows = 1, header = FALSE, sep=sep, stringsAsFactors = FALSE)
   first_data <- read.table(filepath, nrows=1, sep=sep, skip=1)
   if(length(header)==length(first_data)){
      cols <- c("character", rep("numeric", length(header)-1))
   }else if(length(header)==length(first_data)-1){
      cols <- c("character", rep("numeric", length(header)))
   }
   rawdata <- fread(filepath, header=F, sep=sep, skip=1, colClasses=cols)

    ### Again. Add more checking in case there are duplicate rownames
    # read in the first column and check, if there is duplicated rownames (mostly gene names)
    # then send out a warning and then make it unique
    if(sum(duplicated(rawdata$V1))>0){
      warning("There are duplicated rownames in your data")
      warning("Please double check your gene count table")
      rawdata$V1 <- make.names(rawdata$V1, unique=TRUE)
    }

   ### Add more checking in case there are duplicated column names
   # make.names(names, unique=TRUE)
   rawdata %<>% setDF %>% set_rownames(.$V1) %>%
                  '['(colnames(.) != "V1") #%>% as.numeric

   # data doesn't have colnames for first row (rownames)
   if(length(header) == dim(rawdata)[2]){
      # colnames(rawdata) <- unlist(header)
      colnames(rawdata) <- make.names(unlist(header), unique=TRUE)
   }else if (length(header) == dim(rawdata)[2] + 1){
      # colnames(rawdata) <- unlist(header)[-1]
      colnames(rawdata) <- make.names(unlist(header)[-1], unique=TRUE)
   }

   # Add checking data platform
   if(check.platform){
      rawdata <- detect_genecnt_platform(rawdata)
   }

   return(rawdata)
}


rpm <- function(data, mapped_reads){
  data <- t(t(data)/mapped_reads*1000000)
}

run_DESeq2 <- function(countData, annot, column="Time"){
   require(DESeq2)
   register(BiocParallel::MulticoreParam(8))
   group <- dplyr::pull(annot, column) %>% factor()
   coldat <- DataFrame(grp = group)
   ddsfeatureCounts <- DESeqDataSetFromMatrix(countData = countData,
                                              colData = coldat, 
                                              design = ~ grp)
   dds <- DESeq2::DESeq(ddsfeatureCounts, parallel=T)
   detach("package:DESeq2")
   return(dds)
}
  
run_DESeq2_norep <- function(countData, annot, column="Time"){
   require(DESeq2)
   register(BiocParallel::MulticoreParam(8))
   group <- dplyr::pull(annot, column) %>% factor()
   coldat <- DataFrame(grp = group)
   ddsfeatureCounts <- DESeqDataSetFromMatrix(countData = countData,
                                              colData = coldat, 
                                              design = ~ grp)
   rld <- rlogTransformation( dds )
   res <- data.frame(assay(rld), 
                     avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
                     rLogFC = assay(rld)[,2] - assay(rld)[,1] )
   #dds <- DESeq2::DESeq(ddsfeatureCounts, parallel=T)
   detach("package:DESeq2")
   return(res)
}

run_DESeq <- function(countData, annot, column="Time"){
   require(DESeq)
   group <- dplyr::pull(annot, column) %>% factor()
   y <- newCountDataSet(countData, group)
   y <- estimateSizeFactors(y)
   y <- estimateDispersions(y, fitType="local")
   # result <- nbinomTest(y, 2, 1)
   detach("package:DESeq")
   return(y)
}

run_DESeq_norep <- function(countData, annot, column="Time"){
   require(DESeq)
   group <- dplyr::pull(annot, column) %>% factor()
   y <- newCountDataSet(countData, group)
   y <- estimateSizeFactors(y)
   y <- estimateDispersions(y, fitType="local", method="blind",sharingMode="fit-only")
   # result <- nbinomTest(y, 2, 1)
   detach("package:DESeq")
   return(y)
}

run_edgeR <- function(countData, annot, column="Time", method="TMM"){
  require(edgeR)
  group <- dplyr::pull(annot, column) %>% factor()
  y <- DGEList(counts=countData, 
               group=group,
               genes=rownames(countData))
  y <- calcNormFactors(y)
  y <- estimateDisp(y, method=method)
  #y <- estimateTagwiseDisp(y)
  detach("package:edgeR")
  return(y)
}

run_edgeR_norep <- function(countData, annot, column="Time", method="TMM"){
  require(edgeR)
  group <- dplyr::pull(annot, column) %>% factor()
  y <- DGEList(counts=countData, 
               group=group,
               genes=rownames(countData))
  y <- calcNormFactors(y)
  y <- estimateDisp(y, method=method)
  #y <- estimateTagwiseDisp(y)
  detach("package:edgeR")
  return(y)
}

run_voom <- function(countData, annot, column="Time"){
   require(edgeR)
   y <- DGEList(counts=countData)
   y <- calcNormFactors(y) # Scale normalization
   #group <- dplyr::pull(annot, column) %>% rev() %>% factor() 
   group <- dplyr::pull(annot, column) %>% factor(., levels=unique(.) %>% rev)
   design <- model.matrix(~ group) 
   v <- voom(y, design, plot=TRUE) 
   fit <- lmFit(v,design)
   fit <- eBayes(fit)
   detach("package:edgeR")
   return(fit)
}

run_voom_norep <- function(countData, annot, column="Time"){
   require(edgeR)
   y <- DGEList(counts=countData)
   y <- calcNormFactors(y) # Scale normalization
   group <- dplyr::pull(annot, column) %>% rev() %>% factor() 
   #design <- model.matrix(~ group) 
   v <- voom(y, plot=TRUE) 
   fit <- lmFit(v)
   fit <- eBayes(fit)
   detach("package:edgeR")
   return(fit)
}

run_DESeq2_sizefctr <- function(countData, annot, column="Time", control.idx=NULL){
   require(DESeq2)
   register(BiocParallel::MulticoreParam(8))
   group <- dplyr::pull(annot, column) %>% factor()
   coldat <- DataFrame(grp = group)
   ddsfeatureCounts <- DESeqDataSetFromMatrix(countData = countData,
                                              colData = coldat, 
                                              design = ~ grp)
   if(is.numeric(control.idx)){
      if(sum(is.na(countData[control.idx, ])) == 0){
         dds <- estimateSizeFactors(ddsfeatureCounts, controlGenes=house.idx)
         dds <- estimateDispersions(dds)
         dds <- nbinomWaldTest(dds) 
      }else{
         warnings("### Not all the index in control.idx are found")
      }
   }else{
      dds <- DESeq2::DESeq(ddsfeatureCounts, parallel=T)
   }
   detach("package:DESeq2")
   return(dds)
}

run_edgeR_sizefctr <- function(countData, annot, column="Time", method="TMM"){
  require(edgeR)
  group <- dplyr::pull(annot, column) %>% factor()
  y <- DGEList(counts=countData, 
               group=group,
               genes=rownames(countData))
  y <- calcNormFactors(y, method=method)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  detach("package:edgeR")
  return(y)
}


# summary for volcane, logfc, and waterfall plots
sgRNA_summary_plot <- function(deseq2.res, prefix = "sgRNA", n.water.genes=30){
  deseq2.plot.res <- deseq2.res %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(Gene = rownames(deseq2.res), 
           log2FC = log2FoldChange) %>% 
    dplyr::select(Gene, log2FC, padj)
  
  pdf(paste0(prefix, ".Volcano.pdf"), height=10, width=10)
  plot_de_res_volcano(deseq2.plot.res, title = prefix)
  dev.off()
  
  pdf(paste0(prefix, ".log2FC.pdf"), height=10, width=10, title=prefix)
  plot_de_res_log2FC(deseq2.plot.res, ymin = -12, ymax=10)
  dev.off()
  
  tmp <- deseq2.plot.res %>% 
     filter(padj <= 0.5) %>% 
     arrange(log2FC) 
  if(nrow(tmp) >= 30){
     pdf(paste0(prefix, ".waterfall.pdf"), height=8, width=10, title=prefix)
     genelist <- bind_rows(head(tmp, n.water.genes), tail(tmp, n.water.genes)) %>% pull(Gene)
     plot_de_res_waterfall(deseq2.plot.res, genelist=genelist)
     dev.off()
  }else{
     warnings("### Not enough genes pass through the p.adj filter")
  }
}

# volcano
plot_de_res_volcano <- function(de.res, n.genes=10, color="red", title=""){
  # select top n padj, top10 log2FC, top10 -log2FC genes to label
  filter.res <- bind_rows(top_n(de.res, n.genes, -padj),
                          top_n(filter(de.res, is.finite(log2FC)), n.genes, log2FC),
                          top_n(filter(de.res, is.finite(log2FC)), n.genes, -log2FC)) %>% 
    distinct(Gene, .keep_all = TRUE)
  if(nrow(filter.res) >= 100){
     filter.res <- filter.res %>% 
        arrange(abs(log2FC)) %>% 
        dplyr::slice(1:50)
  }
  
  xlim.range <- de.res$log2FC %>% range %>% abs %>% max %>% round() %>% "+"(3)
  r_volcano <- de.res %>% 
    ggplot(aes(x=log2FC, -log10(padj), label=Gene)) + 
    geom_point(alpha=0.6) + 
    geom_point(data = filter.res, col=color, alpha=0.6) + 
    theme_bw() + 
    xlim(-xlim.range, xlim.range) +
    ggtitle(title)
  
  r_volcano_label <- r_volcano + 
    geom_text_repel(data= filter.res, 
                    aes(label=Gene), size=2, col="gray15", 
                    point.padding = 0.6) 
  
  print(r_volcano)
  print(r_volcano_label)
}


plot_de_res_log2FC <- function(de.res, n.genes=10, color="red", title="", ymin=-6, ymax=6){
  # color <- "red"
  # n.genes <- 10
  de.res <- arrange(de.res, log2FC) %>% 
    mutate(index = 1:nrow(.)) 
  
  sig.res <- top_n(de.res, n.genes, -padj)
  if(nrow(sig.res) >= 30){
    sig.res <- sig.res %>% 
       arrange(abs(log2FC)) %>% 
       dplyr::slice(1:30)
  }
  
  top.res <- bind_rows(top_n(filter(de.res, is.finite(log2FC)), n.genes, log2FC),
                       top_n(filter(de.res, is.finite(log2FC)), -n.genes, log2FC)) %>% 
    distinct(Gene, .keep_all = TRUE)
  if(nrow(top.res) >= 30){
    top.res <- top.res %>% 
       arrange(abs(log2FC)) %>% 
       dplyr::slice(1:30)
  }
  
  r_log2FC <- de.res %>% 
    ggplot(aes(x=index, y=log2FC)) + 
    geom_point(size=1, alpha=0.5) + 
    geom_point(data = top.res, col="orange", alpha=0.6) + 
    geom_point(data = sig.res, col=color, alpha=0.9) + 
    theme_classic() + 
    geom_hline(yintercept = c(0), linetype="solid", lwd=0.2) +
    geom_hline(yintercept = c(-1, 1), linetype="longdash", lwd=0.2) + 
    ylim(c(ymin, ymax)) +
    xlab("") +
    ylab(paste0("log2FoldChange( ", prefix, " )")) +
    ggtitle(title)
  r_log2FC_label <- r_log2FC + 
    geom_text_repel(data= sig.res, 
                    aes(label=Gene), size=3, col=color,
                    point.padding = 1) + 
    geom_text_repel(data= top.res, 
                    aes(label=Gene), size=2, col="gray15", 
                    point.padding = 1) 
  
  print(r_log2FC)
  print(r_log2FC_label)
}
  
  
# random pick 100 genes
plot_de_res_waterfall <- function(de.res, genelist, title="", size=3){
   water.de.res <- filter(de.res, Gene %in% genelist) %>% 
    arrange(log2FC) %>% 
    mutate(index = 1:nrow(.)) %>% 
    #mutate(Group = gsub("-.*", "", Gene)) 
    mutate(Group = gsub("\\..*", "", Gene)) 
  
  ylim.range <- water.de.res$log2FC %>% abs %>% max %>% round() %>% "+"(3)
  
  r_waterfall <- water.de.res %>% 
    ggplot(aes(x=index, y=log2FC)) + 
    geom_bar(stat="identity", 
             width=0.7, 
             position = position_dodge(width=0.4)) +
    ggtitle(title) + 
    xlab("") + 
    ylab(paste0("log2FoldChange( ", prefix, " )")) +
    theme_classic() + 
    theme(axis.line.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=size + 10),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(face="bold",angle=90)) +
    coord_cartesian(ylim = c(-ylim.range, ylim.range))
  r_waterfall_label1 <- r_waterfall + 
    geom_text(data = subset(water.de.res, log2FC >= 0), aes(label=Gene), angle = 90, vjust=0.5, hjust=-0.2, size=size) + 
    geom_text(data = subset(water.de.res, log2FC < 0), aes(label=Gene), angle = 90, vjust=0.5, hjust=1.2, size=size) 
  
  r_waterfall <- water.de.res %>% 
    ggplot(aes(x=index, y=log2FC, fill=Group)) + 
    geom_bar(stat="identity", 
             width=0.7, 
             position = position_dodge(width=0.4)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(length(genelist))) + 
    ggtitle(title) + 
    xlab("") + 
    ylab(paste0("log2FoldChange( ", prefix, " )")) +
    theme_classic() + 
    theme(axis.line.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=size + 10),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(face="bold",angle=90)) +
    theme(legend.position = "none") + 
    coord_cartesian(ylim = c(-ylim.range, ylim.range))
  r_waterfall_label2 <- r_waterfall + 
    geom_text(data = subset(water.de.res, log2FC >= 0), aes(label=Gene), angle = 90, vjust=0.5, hjust=-0.2, size=size) + 
    geom_text(data = subset(water.de.res, log2FC < 0), aes(label=Gene), angle = 90, vjust=0.5, hjust=1.2, size=size) 

  print(r_waterfall_label1) 
  print(r_waterfall_label2) 
}




convert_tibble_to_dataframe <- function(x){
  res <- as.data.frame(x[, -1])
  rownames(res) <- x[, 1] %>% unlist()
  return(res)
}

convert_data_to_tibble <- function(x, name="Gene"){
  res <- as.data.frame(x) %>% 
     rownames_to_column(name) %>% 
     as_tibble()
  return(res)
}

pheatmap_fontsize_row <- function(x){
   if(class(x) %in% c("data.frame", "matrix")){
      n.genes <- nrow(x)
   }else if(class(x) == "character"){
      n.genes <- length(x)
   }else{
      stop(paste("unknown data type", class(x)))
   }

   fontsize_row <- case_when(
                             n.genes > 1000 ~ 0.5,
                             n.genes > 800 ~ 0.6,
                             n.genes > 500 ~ 1,
                             n.genes > 300 ~ 2,
                             n.genes > 200 ~ 3,
                             n.genes > 120 ~ 4,
                             n.genes > 90 ~ 5,
                             n.genes > 60 ~ 6,
                             n.genes > 40 ~ 7,
                             n.genes > 20 ~ 8,
                             n.genes > 10 ~ 9,
                             TRUE ~ 10
                             )
   return(fontsize_row)
}

pheatmap_fontsize_col <- function(x){
   if(class(x) %in% c("data.frame", "matrix")){
      n.samples <- ncol(x)
   }else if(class(x) == "character"){
      n.samples <- length(x)
   }else{
      stop(paste("unknown data type", class(x)))
   }

   fontsize_col <- case_when(
                             n.samples > 1000 ~ 0.4,
                             n.samples > 800 ~ 0.6,
                             n.samples > 500 ~ 1,
                             n.samples > 300 ~ 2,
                             n.samples > 200 ~ 3,
                             n.samples > 120 ~ 4,
                             n.samples > 100 ~ 4.5,
                             n.samples > 80 ~ 5,
                             n.samples > 60 ~ 5.5,
                             n.samples > 40 ~ 6,
                             n.samples > 20 ~ 7,
                             TRUE ~ 8
                             )
   return(fontsize_col)
}

pheatmap_pdf_width <- function(x){
   if(class(x) %in% c("data.frame", "matrix")){
      n.samples <- ncol(x)
   }else if(class(x) == "character"){
      n.samples <- length(x)
   }else{
      stop(paste("unknown data type", class(x)))
   }

   pdf_width <- case_when(
                          n.samples > 500 ~ 22,
                          n.samples > 400 ~ 20,
                          n.samples > 350 ~ 19,
                          n.samples > 300 ~ 18,
                          n.samples > 250 ~ 17,
                          n.samples > 200 ~ 16,
                          n.samples > 150 ~ 15,
                          n.samples > 120 ~ 14,
                          n.samples > 90 ~ 13,
                          n.samples > 70 ~ 12,
                          n.samples > 50 ~ 11,
                          n.samples > 30 ~ 10,
                          n.samples > 24 ~ 9,
                          n.samples > 17 ~ 8,
                          n.samples > 10 ~ 7.5,
                          TRUE ~ 7
                          )
   return(pdf_width)
}

pheatmap_pdf_height <- function(x){
   if(class(x) %in% c("data.frame", "matrix")){
      n.genes <- nrow(x)
   }else if(class(x) == "character"){
      n.genes <- length(x)
   }else{
      stop(paste("unknown data type", class(x)))
   }

   pdf_height <- case_when(
                           n.genes > 300 ~ 20,
                           n.genes > 100 ~ 15,
                           n.genes > 80 ~ 12,
                           n.genes > 70 ~ 11,
                           n.genes > 60 ~ 10,
                           n.genes > 50 ~ 9,
                           n.genes > 40 ~ 10,
                           n.genes > 30 ~ 9,
                           n.genes > 25 ~ 8.5,
                           n.genes > 20 ~ 8,
                           n.genes > 15 ~ 7.5,
                           n.genes > 12 ~ 7,
                           n.genes > 9 ~ 6.5,
                           n.genes > 7 ~ 6,
                           n.genes > 5 ~ 5,
                           n.genes > 3 ~ 4.5,
                           n.genes > 2 ~ 4.5,
                           TRUE ~ 4 
                           )
   return(pdf_height)
}

# divide always works row-wise
# first time to divide each transcript by length (row-wise)
# transpose second time to divide samples by sample_total_reads 
cnt_to_fpkm <- function(x, gene_len){
   num.reads <- colSums(x)
   res <- t(t(x/gene_len)/num.reads) * 10^9
   return(res)
}

# transcript first normalized by length
# sum the total normalized reads
# calculate the proportionn as tpm
cnt_to_tpm <- function(x, gene_len){
   x_norm_len <- x/gene_len
   sum_norm_len <- colSums(x_norm_len)
   res <- t(t(x_norm_len)/sum_norm_len) * 10^6
}
  




stat_smooth_func <- function(mapping = NULL, data = NULL,
                        geom = "smooth", position = "identity",
                        ...,
                        method = "auto",
                        formula = y ~ x,
                        se = TRUE,
                        n = 80,
                        span = 0.75,
                        fullrange = FALSE,
                        level = 0.95,
                        method.args = list(),
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE,
                        xpos = NULL,
                        ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,
                      
                      setup_params = function(data, params) {
                        # Figure out what type of smoothing to do: loess for small datasets,
                        # gam with a cubic regression basis for large data
                        # This is based on the size of the _largest_ group.
                        if (identical(params$method, "auto")) {
                          max_group <- max(table(data$group))
                          
                          if (max_group < 1000) {
                            params$method <- "loess"
                          } else {
                            params$method <- "gam"
                            params$formula <- y ~ s(x, bs = "cs")
                          }
                        }
                        if (identical(params$method, "gam")) {
                          params$method <- mgcv::gam
                        }
                        
                        params
                      },
                      
                      compute_group = function(data, scales, method = "auto", formula = y~x,
                                               se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                               xseq = NULL, level = 0.95, method.args = list(),
                                               na.rm = FALSE, xpos=NULL, ypos=NULL) {
                        if (length(unique(data$x)) < 2) {
                          # Not enough data to perform fit
                          return(data.frame())
                        }
                        
                        if (is.null(data$weight)) data$weight <- 1
                        
                        if (is.null(xseq)) {
                          if (is.integer(data$x)) {
                            if (fullrange) {
                              xseq <- scales$x$dimension()
                            } else {
                              xseq <- sort(unique(data$x))
                            }
                          } else {
                            if (fullrange) {
                              range <- scales$x$dimension()
                            } else {
                              range <- range(data$x, na.rm = TRUE)
                            }
                            xseq <- seq(range[1], range[2], length.out = n)
                          }
                        }
                        # Special case span because it's the most commonly used model argument
                        if (identical(method, "loess")) {
                          method.args$span <- span
                        }
                        
                        if (is.character(method)) method <- match.fun(method)
                        
                        base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                        model <- do.call(method, c(base.args, method.args))
                        
                        m = model
                        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                         list(a = format(coef(m)[1], digits = 3), 
                                              b = format(coef(m)[2], digits = 3), 
                                              r2 = format(summary(m)$r.squared, digits = 3)))
                        func_string = as.character(as.expression(eq))
                        
                        if(is.null(xpos)) xpos = min(data$x)*0.9
                        if(is.null(ypos)) ypos = max(data$y)*0.9
                        data.frame(x=xpos, y=ypos, label=func_string)
                        
                      },
                      
                      required_aes = c("x", "y")
)


lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
      b = format(abs(coef(m)[2]), digits = 2),
      r2 = format(summary(m)$r.squared, digits = 3));
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  as.character(as.expression(eq));                 
}


annotate_data <- function(rawdata, gene_col="Gene_name", species="Human"){
  if(species == "Human"){
     annot_file <- "~/Copy/MSKCC/gtf_files/GRCh38_Ensembl91_p10_mart_export_GO_20180308.txt"
  }else if(species == "Mouse"){
     annot_file <- "~/Copy/MSKCC/gtf_files/GRCm38_Ensembl91_p5_mart_export_GO_20180308.txt"
  }else{
     warning("### Error: we only support 'human' and 'mouse' for now")
     exit
  }

  # for each of the ensembl_gene_id, there are multiple transcripts
  # we will summarise all the transcript_length and use that for fpkm/tpm calculation
  GRC38.91.ensembl <- read_tsv(annot_file) %>%
    set_names(gsub(" ", "_", colnames(.))) %>% 
    rename(transcript_length = `Transcript_length_(including_UTRs_and_CDS)`,
           chromosome_name = `Chromosome/scaffold_name`, 
           start_position = `Gene_start_(bp)`,
           end_position = `Gene_end_(bp)`
    ) %>% 
    dplyr::select(-contains("Transcript", ignore.case = FALSE))
  
  # we have to do it for transcript_length and GO_term separately
  GRC38.91.ensembl.length <- GRC38.91.ensembl %>% 
    dplyr::select(Gene_stable_ID, transcript_length) %>%
    group_by(Gene_stable_ID) %>% 
    summarise_all(sum) #%>% 
  GRC38.91.ensembl.GO <- GRC38.91.ensembl %>% 
    dplyr::select(Gene_stable_ID, GO_term_name, GO_term_accession) %>% 
    group_by(Gene_stable_ID) %>% 
    summarise(GO_term = paste(ifelse(is.na(GO_term_name), "", GO_term_name), collapse = "/"),
              GO_accession = paste(ifelse(is.na(GO_term_accession), "", GO_term_accession), collapse = "/")) %>% 
    mutate(GO_term = gsub("^\\/{1,}$", NA, GO_term)) %>% 
    mutate(GO_term = gsub("^\\/{1,}", "", GO_term)) %>% 
    mutate(GO_term = gsub("\\/{1,}$", "", GO_term)) %>% 
    mutate(GO_accession = gsub("^\\/{1,}$", NA, GO_accession)) %>% 
    mutate(GO_accession = gsub("^\\/{1,}", "", GO_accession)) %>% 
    mutate(GO_accession = gsub("\\/{1,}$", "", GO_accession)) %>% 
    mutate(GO_term = ifelse(nchar(GO_term) == 0, NA, GO_term), 
           GO_accession = ifelse(nchar(GO_accession) == 0, NA, GO_accession))
  
  GRC38.91.ensembl <-
    left_join(dplyr::select(GRC38.91.ensembl, -transcript_length, -contains("GO")) %>% distinct(.), GRC38.91.ensembl.length) %>% 
    left_join(., GRC38.91.ensembl.GO, by="Gene_stable_ID") %>%
    dplyr::select(Gene_stable_ID, Gene_name:end_position, transcript_length, everything()) %>%
    mutate(Strand = ifelse(Strand == 1, "+", "-")) %>% 
    arrange(Gene_stable_ID)
  
  # match the order of the info
  # res <- left_join(rawdata, GRC38.91.ensembl, by=c( (!!rlang::sym(gene_col))="Gene_stable_ID"))
   res <- left_join(rawdata, GRC38.91.ensembl, by="Gene_name")
  #res <- left_join(rawdata, GRC38.91.ensembl, by=c("ID"="Gene_name"))
  return(res)
}


annotate_data_HPA_protein <- function(rawdata, gene_col="Gene", species="Human", simple=FALSE){
   protein <- read_tsv("~/Copy/MSKCC/Human_Protein_Atlas/proteinatlas.tsv") %>% 
      set_names(gsub(" ", "_", colnames(.))) 
   if(simple){
      protein <- protein %>% 
         dplyr::select(Gene, Gene_description:Antibody, Subcellular_location, `TPM_max_in_non-specific`) 
   }
   # match the order of the info
   res <- left_join(rawdata, protein, by="Gene")

   return(res)
}

# add prognostic info
annotate_data_HPA_pathology <- function(rawdata, tissue="pancreatic"){
   pathology <- read_tsv("~/Copy/MSKCC/Human_Protein_Atlas/pathology.tsv") %>% 
      set_names(gsub(" - ", "_", colnames(.))) %>% 
      set_names(gsub(" ", "_", colnames(.))) %>% 
      mutate(Cancer = gsub(" cancer", "", Cancer)) %>% 
      filter(grepl(tissue, Cancer)) %>% 
      dplyr::select(Gene= Gene_name, Anti_High = High, Anti_Medium = Medium, Anti_Low = Low, 
                    Anti_Not_detected = Not_detected, prognostic_favourable, prognostic_unfavourable)

   if(nrow(pathology) > 0){
      res <- left_join(rawdata, pathology, by="Gene")
      return(res)
   }else{
      print(paste("Something wrong with the tissue name:", tissue, "please check again"))
      return(NA)
   }
}

annotate_data_HPA_tissue <- function(rawdata, tissue="pancreas"){
   # add tissue specific TPM
   tis <- read_tsv("~/Copy/MSKCC/Human_Protein_Atlas/rna_tissue.tsv") %>% 
      set_names(gsub(" ", "_", colnames(.))) %>% 
      filter(Sample == tissue) %>% 
      dplyr::select(Gene = Gene_name, Value) %>% 
      dplyr::rename(!!sym(paste0(tissue, "_TPM")) := Value)

   if(nrow(tis) > 0){
      res <- left_join(rawdata, tis, by="Gene")
      return(res)
   }else{
      print(paste("Something wrong with the tissue name:", tissue, "; please check again"))
      return(NA)
   }
}

annotate_data_addlink <- function(rawdata, gene_col="Gene"){
   link <- read_csv("~/Copy/MSKCC/gtf_files/GRCh38_Ensembl91_p10_mart_export_AnnotLinks.csv") %>% 
      dplyr::select(-Gene_stable_ID) %>% 
      group_by(Gene) %>% 
      dplyr::slice(1) %>% 
      ungroup()
   res <- left_join(rawdata, link, by="Gene")
   return(res)
}


annotate_data_by_ensemblID <- function(rawdata, id_col="Gene_stable_ID", species="Human"){
  if(species == "Human" | species == "human"){
     annot_file <- "~/Copy/MSKCC/gtf_files/GRCh38_Ensembl91_p10_mart_export_GO_20180308.txt"
  }else if(species == "Mouse" | species == "mouse"){
     annot_file <- "~/Copy/MSKCC/gtf_files/GRCm38_Ensembl91_p5_mart_export_GO_20180308.txt"
  }else{
     warning("### Error: we only support 'human' and 'mouse' for now")
     exit()
  }

  # for each of the ensembl_gene_id, there are multiple transcripts
  # we will summarise all the transcript_length and use that for fpkm/tpm calculation
  GRC38.91.ensembl <- read_tsv(annot_file) %>%
    set_names(gsub(" ", "_", colnames(.))) %>% 
    rename(transcript_length = `Transcript_length_(including_UTRs_and_CDS)`,
           chromosome_name = `Chromosome/scaffold_name`, 
           start_position = `Gene_start_(bp)`,
           end_position = `Gene_end_(bp)`
    ) %>% 
    dplyr::select(-contains("Transcript", ignore.case = FALSE))
  
  # we have to do it for transcript_length and GO_term separately
  GRC38.91.ensembl.length <- GRC38.91.ensembl %>% 
    dplyr::select(Gene_stable_ID, transcript_length) %>%
    group_by(Gene_stable_ID) %>% 
    summarise_all(sum) #%>% 
  GRC38.91.ensembl.GO <- GRC38.91.ensembl %>% 
    dplyr::select(Gene_stable_ID, GO_term_name, GO_term_accession) %>% 
    group_by(Gene_stable_ID) %>% 
    summarise(GO_term = paste(ifelse(is.na(GO_term_name), "", GO_term_name), collapse = "/"),
              GO_accession = paste(ifelse(is.na(GO_term_accession), "", GO_term_accession), collapse = "/")) %>% 
    mutate(GO_term = gsub("^\\/{1,}$", NA, GO_term)) %>% 
    mutate(GO_term = gsub("^\\/{1,}", "", GO_term)) %>% 
    mutate(GO_term = gsub("\\/{1,}$", "", GO_term)) %>% 
    mutate(GO_accession = gsub("^\\/{1,}$", NA, GO_accession)) %>% 
    mutate(GO_accession = gsub("^\\/{1,}", "", GO_accession)) %>% 
    mutate(GO_accession = gsub("\\/{1,}$", "", GO_accession)) %>% 
    mutate(GO_term = ifelse(nchar(GO_term) == 0, NA, GO_term), 
           GO_accession = ifelse(nchar(GO_accession) == 0, NA, GO_accession))
  
  GRC38.91.ensembl <-
    left_join(dplyr::select(GRC38.91.ensembl, -transcript_length, -contains("GO")) %>% distinct(.), GRC38.91.ensembl.length) %>% 
    left_join(., GRC38.91.ensembl.GO, by="Gene_stable_ID") %>%
    dplyr::select(Gene_stable_ID, Gene_name:end_position, transcript_length, everything()) %>%
    mutate(Strand = ifelse(Strand == 1, "+", "-")) %>% 
    arrange(Gene_stable_ID)
  
  # match the order of the info
  # res <- left_join(rawdata, GRC38.91.ensembl, by=c( (!!rlang::sym(gene_col))="Gene_stable_ID"))
   res <- left_join(rawdata, GRC38.91.ensembl, by="Gene_stable_ID")
  #res <- left_join(rawdata, GRC38.91.ensembl, by=c("ID"="Gene_name"))
  return(res)
}


donuts <- function(x, group = 1, labels = NA, col = NULL, radius = c(.7, 1)) {
  group <- rep_len(group, length(x))
  ug  <- unique(group)
  tbl <- table(group)[order(ug)]

  col <- if (is.null(col))
    seq_along(ug) else rep_len(col, length(ug))
  col.main <- Map(rep, col[seq_along(tbl)], tbl)
  col.sub  <- lapply(col.main, function(x) {
    al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    Vectorize(adjustcolor)(x, alpha.f = al)
  })

  plot.new()

  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L],
      col = unlist(col.sub), labels = labels)

  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L],
      col = unlist(col.main), labels = NA)
}


mydonuts <- function(x, group = 1, labels = NA, col = NULL, radius = c(.7, 1), clockwise=TRUE, main=NULL) {
  group <- rep_len(group, length(x))
  ug  <- unique(group)
  tbl <- table(group)[order(ug)]

  col <- if (is.null(col))
    seq_along(ug) else rep_len(col, length(ug))
  col.main <- Map(rep, col[seq_along(tbl)], tbl)
  col.sub  <- lapply(col.main, function(x) {
    al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    Vectorize(adjustcolor)(x, alpha.f = al)
  })

  plot.new()

  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L], 
      clockwise = clockwise, main = main,
      col = unlist(col.sub), labels = labels)

  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L], clockwise = clockwise,
      col = unlist(col.main), labels = NA)
}


parse_file_as_list <- function(dat){
   nams <- dat[, 1]
   dat <- dat[, -1] 
   ldat <- split(dat, seq_len(nrow(dat)))
   ldat <- lapply(ldat, function(x) x[!is.na(x)][-1]) # remove the first http weblink
   names(ldat) <- nams

   return(ldat)
}

count_upper_rows <- function(x){
   num.upper <- apply(data.frame(x), 1, function(x) gsub("\\d", "", x)) %>% 
      as.data.frame() %>% 
      apply(., 1, function(x) gsub("[-.]", "", x)) %>% 
      as.data.frame() %>% 
      apply(., 1, function(x) gsub("orf", "", x)) %>% 
      as.data.frame() %>% 
      apply(., 1, function(x) grepl("^[[:upper:]]+$", x)) %>% 
      sum
   return(num.upper)
}


check_genes_species <- function(x){
}

convert_human_to_mouse <- function(gene.list, colname="Gene"){
   human.gene <- read_csv("~/Copy/MSKCC/gtf_files/Ensembl91_mouse_to_human_gene.csv")

   res <- dplyr::select(gene.list, !!sym(colname)) %>% 
      left_join(., human.gene, by=c("Gene"="HGNC.symbol")) %>% 
      distinct(Gene, .keep_all = TRUE) %>% 
      left_join(gene.list, ., by="Gene") %>% 
      mutate(MGI.symbol = ifelse(is.na(MGI.symbol), !!sym(colname), MGI.symbol)) %>% 
      dplyr::select(-Gene) %>%  
      dplyr::select(Gene = MGI.symbol, everything()) 
   return(res)
}

convert_mouse_to_human <- function(gene.list, colname="Gene"){
   human.gene <- read_csv("~/Copy/MSKCC/gtf_files/Ensembl91_mouse_to_human_gene.csv")

   res <- dplyr::select(gene.list, !!sym(colname)) %>% 
      left_join(., human.gene, by=c("Gene"="MGI.symbol")) %>% 
      distinct(Gene, .keep_all = TRUE) %>% 
      left_join(gene.list, ., by="Gene") %>% 
      mutate(HGNC.symbol = ifelse(is.na(HGNC.symbol), !!sym(colname), HGNC.symbol)) %>% 
      dplyr::select(-Gene) %>%  
      dplyr::select(Gene = HGNC.symbol, everything())
   return(res)
}


remove_constant_row <- function(dat){
   idx <- apply(dat, 1, range) %>% 
      apply(., 2, function(x) x[1] == x[2]) %>% 
      which
   if(length(idx) > 0){
      dat <- dat[-idx, ]
   }
   return(dat)
}


#' Check rownames from data and auto detect which platform it is from
#'
#' Basically we want to convert those affy metrix array data and just
#' leave normal gene_symbols intact
#' @param data Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' detect_genecnt_platform("TP53")

detect_genecnt_platform <- function(data, method="mean"){
	require("annotate")
	require(magrittr)

	PROBES <- as.character(rownames(data))
	num.genes <- nrow(data)

	# for now just use simple nrow as check point
	if(num.genes == 22283){
  	  require(hgu133a.db) # 22283 probes
	  GPL="GPL96"
	  print ("Guess it's hgu133a")
	}else if(num.genes == 22645){
  	  require(hgu133b.db) # 22645 probes
	  GPL="GPL97"
	  print ("Guess it's hgu133b")
	}else if(num.genes == 54675 || num.genes == 44137){
  	  require(hgu133plus2.db) # 54675 probes
  	  # hard fix here for now, need to add mapped rownames to all these probe.annotation and then decide which one 
	  GPL="GPL570"
	  print ("Guess it's hgu133plus2")
	}else if(num.genes == 22277){
  	  require(hgu133a2.db) 
	  GPL="GPL571"
	  print ("Guess it's hgu133a2")
	}else{
	  print ("Guess its normal RNA-Seq, do nothing")
	  return(data)
	}

	if(method == "mean"){
	  select.fun <- function(x) rowMeans(x)
	}else if(method == "sd"){
	  select.fun <- function(x) rowSds(x)
	}else if(method == "max"){
	  select.fun <- function(x) rowMax(x)
	}else if(method == "min"){
	  select.fun <- function(x) rowMin(x)
	}else if(method == "iqr"){
	  select.fun <- function(x) apply(x, 1, IQR)
	}
	# There are two levels of mutliple mapping
	# one is that each probe_ids may have multiple mapped gene_symbols
	# For now, we just select the first gene_symbols if there are multiples
	gene.symbols <- affy_probe_to_gene_symbol(PROBES, GPL) %>%
	                filter(!duplicated(PROBEID)) %>%
	                dplyr::select(SYMBOL) %>%
	                unlist

	# Another place is that there may be multiple probes designed for each gene
	# Now we select the probes with the max(mean expression across samples)
	# Mayb fix it later by adding more different filter selections? Ray. 2015-10-30
	data.gene <- data %>% as.data.table %>%
	                      '['(, RowExp := select.fun(.SD)) %>%
	                      '['(, SYMBOL := gene.symbols) %>% group_by(., SYMBOL) %>%
	                      filter(., which.max(RowExp))

	# Remove rows without gene_symbols (NAs), then convert back to data.frame
	final.data <- data.gene %>%
	                  filter(!is.na(SYMBOL)) %>%
	                  setDF %>%
	                  set_rownames(.$SYMBOL) %>%
	                  '['(, !(colnames(.) %in% c("RowExp", "SYMBOL")))

	return(final.data)
}

#' Check rownames from data and auto detect which platform it is from
#'
#' Basically we want to convert those affy metrix array data and just
#' leave normal gene_symbols intact
#' @param data Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' detect_genecnt_platform("TP53")

detect_genecnt_ID <- function(data, sum.method="mean"){
	require(magrittr)
	require(pathview)

	# Map molecular data onto Entrez Gene IDs
	id.map.sym2eg <- id2eg(ids = rownames(data), category = "SYMBOL", org="Hs")
	data.entrez <- mol.sum(mol.data = data, id.map = id.map.sym2eg, sum.method = sum.method)

	# deseq2.fc <- gene.entrez[, 1]
	# exp.fc=deseq2.fc
	# head(exp.fc)
	# out.suffix="deseq2"
	return(data.entrez)
}

#' This function allows you to run DESeq normalization on your data 
#' @param countdata integer data frame
#' @keywords Normalization
#' @export
#' @examples
#' DESeq_Normalization(countdata)
DESeq_Normalization <- function(countdata){
   	require(DESeq)
   	condition <- factor(rep("Sample", ncol(countdata)))
   	countdata <- newCountDataSet(countdata,condition )
   	countdata <- estimateSizeFactors( countdata )
   	normalized.countdata <- counts(countdata, normalized=TRUE)
   	#return(exprs(normalized.countdata))
   	return(normalized.countdata)
}

#' Extract data by MAD value 
#'
#' This function allows you to express your love of cats.
#' @param data Original gene expression data 
#' @param topN How many genes from TopMAD list
#' @keywords cats
#' @export
#' @examples
#' cat_function()
extract_data_by_mad <- function (data, topN, by="row", type="data"){
	by.idx <- ifelse(by=="row", 1, 2)
	data.mad.genes <- apply(data, by.idx, mad) %>% 
					  select.top.n(., topN, bottom=F) %>% 
					  names
	if(type=="data"){
		if(by=="row"){
			data <- data[rownames(data) %in% data.mad.genes, ]
		}else if(by=="col"){
			data <- data[ ,colnames(data) %in% data.mad.genes]
		}else{
			stop("Wrong parameters, by can be 'row' or 'col'")
		}
		return(data)
	}else if(type=="genes"){
		return(data.mad.genes)
	}else{
		stop("Wrong parameters, type can be 'data' or 'genes'")
	}
}

#' Select top number from a vector
#'
#' Select top # from a vector
#' @param data numerical data
#' @param n top number
#' @param bottom whether to select from bottom?
#' @keywords select 
#' @export
#' @examples
#' select.top.n(data, 100, bottom=F)

select.top.n <- function(data, n, bottom=F){
   if(bottom){
      data <- head(sort(data, decreasing = F), n=n)
   }else{
      data <- head(sort(data, decreasing = T), n=n)
   }
   return(data)
}


#' Convert name to color
#'
#' This function convert the sample_name in the data into different colors 
#' @param name usually sample name, separate with "_" 
#' @keywords color
#' @export
#' @examples
#' name_to_color()
# convert sample name to color
name_to_color <- function(name, split_pattern="\\_", num_color=1,
						  ColScheme=c("naikai", "naikai2", "Set1", "Set2", "Set3", "Dark2", "Pastel1", "Pastel2", "Paired", "Accent")
						 )
{
	require(magrittr)

	col.color <- list()
	ColSideColors <- NULL
	col.color.name <- list()
	ColSideColors.name <- NULL

	column.names <- strsplit(name, split=split_pattern)
	max_num_breaks <- max(sapply(column.names, length))

	if(num_color != length(ColScheme)){
		stop(paste0("Error: num_color:", num_color, "is different than num of ColScheme:", length(ColScheme), "\n"))
	}
	if(max_num_breaks > num_color){
		warning("There are more possible breaks in the name than the num_color specified. We will use only the num_color")
	}

	for(i in 1:num_color){
		num_names <- sapply(column.names, function(x) x[i]) %>% unique %>% length
		col.color[[i]] <- create.brewer.color(sapply(column.names, function(x) x[i]), num_names, ColScheme[i])
		ColSideColors <- cbind(ColSideColors, col.color[[i]])
		col.color.name[[i]] <- sapply(column.names, function(x) x[i])
		ColSideColors.name <- cbind(ColSideColors.name, col.color.name[[i]])
	}

	res <- list()
	if(num_color==1){
		if(length(unique(as.character(ColSideColors[,1])))==1){
			res[["color"]] <- NULL
			res[["name"]] <- NULL
		}else{
			res[["color"]] <- as.matrix(ColSideColors[, 1])
			res[["name"]] <- as.matrix(ColSideColors.name[, 1])
		}
	}else if(num_color>1){
		num_color <- ifelse(num_color > max_num_breaks, max_num_breaks, num_color)
		res[["color"]] <- as.matrix(ColSideColors[, 1:num_color])
		res[["name"]] <- as.matrix(ColSideColors.name[, 1:num_color])
	}else{
		stop("Error: num_color must be greater than 0")
	}

	return(res)
}


### run enrichR using DE analysis results
run_enrichR <- function(dat, gene_col = "Gene_name", dbs = NULL, do.save = FALSE, min.padj=0.05, prefix = "GO"){
  require(enrichR)
  ### enrichR for sig.DE genes 
  if(is.null(dbs)){
    dbs <- c("GO_Molecular_Function_2017b", 
             "GO_Cellular_Component_2017b", 
             "GO_Biological_Process_2017b",
             "BioCarta_2016",
             "KEGG_2016", 
             "Reactome_2016")
  }
  GO.res <- paste0(prefix, ".rda")
  
  if(file.exists(GO.res)){
     print(paste("### GO result file exists:", GO.res, "skipped running"))
     load(GO.res)
  }else{
     go.res <- tibble(content = list(dat)) %>% 
        dplyr::mutate(both.sig = map(content, ~ extract_sig_genes(.x, dir="both", min.padj=min.padj)),
               up.sig = map(content, ~ extract_sig_genes(.x, dir="up", min.padj=min.padj)),
               down.sig = map(content, ~ extract_sig_genes(.x, dir="down", min.padj=min.padj))) %>%
     dplyr::select(-content) %>%
     dplyr::mutate(both.enrichR = map(both.sig, ~ enrichr(.x, dbs)), 
            up.enrichR = map(up.sig, ~ enrichr(.x, dbs)),
            down.enrichR = map(down.sig, ~ enrichr(.x, dbs))
            )
  }
  
  if(do.save){
    save(go.res, file = GO.res)
  }
  return(go.res)
}



### separate gene list based on it's log2FC
extract_sig_genes <- function(x, dir="both", 
                           gene_col = "Gene_name",
                           log2fc_col = "log2FoldChange", min.log2FC=1, 
                           p_col = "padj", min.padj=0.05){
  if(dir=="both"){
    res <- dplyr::filter(x, abs(!!rlang::sym(log2fc_col)) >= min.log2FC & !!rlang::sym(p_col) <= min.padj) 
  }else if(dir == "up"){
    res <- dplyr::filter(x, !!rlang::sym(log2fc_col) >= min.log2FC & !!rlang::sym(p_col) <= min.padj)
  }else if(dir == "down"){
    res <- dplyr::filter(x, !!rlang::sym(log2fc_col) < min.log2FC & !!rlang::sym(p_col) <= min.padj)
  }else{
    stop("## Don't know what is this dir:", dir, "please try again with 'both', 'up' or 'down'")
  } 

  gene <- res %>% 
   arrange(desc(abs(!!rlang::sym(log2fc_col)))) %>% 
   pull(!!rlang::sym(gene_col))

  return(gene)
}


# extract top genes from enrichR results
extract_top_enrichR_res <- function(dat, top=20, by="Adjusted.P.value"){
  if(by == "Adjusted.P.value"){
    res <- dat %>%
      dplyr::select(Term, score = Adjusted.P.value, Type) %>%
      mutate(score = -log10(score)) 
  }else if(by == "Combined.Score"){
    res <- dat %>%
      dplyr::select(Term, score = Combined.Score, Type)
  }else{
    stop(paste("by:", by, "is not recognised, please try with 'padj' or 'combined.score'"))
  }
  res <- res %>% 
    group_by(Type) %>%
    arrange(desc(score)) %>%
    dplyr::slice(1:top) %>%
    ungroup %>%
    spread(key = Type, value = score) %>%
    replace(is.na(.), 0) %>%
    distinct(Term, .keep_all = TRUE) %>% 
    arrange(desc(!!rlang::sym(colnames(.)[2])))
  
  return(res)
}


enrichR_heatmap <- function(enrichR_top_res, enrich.sig.pdf, height=12, width=12){
   pdf(enrich.sig.pdf, height=height, width=width)

   for(name in names(enrichR_top_res)){
      enrichR_top_res[[name]] %>%
         mutate(Term = gsub("_Homo sapiens.*", "", Term)) %>%
         # mutate(Term = gsub("\\(GO.*", "", Term)) %>%
         convert_tibble_to_dataframe() %>%
         pheatmap(., main = name, color = color, breaks = breaks2, cellwidth=40)
   }
   dev.off()
}

go_barplot <- function(name, dat, save=FALSE, text.size=10, height=8, width=12){
  aa <- dat %>% 
    ggplot(aes(x=Term, y=Score, fill=Pathway)) + 
    theme_classic() + 
    geom_col(position = "dodge")  + 
    coord_flip() +
    theme(legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_fill_brewer(palette = "Set1") + 
    ggtitle(name) + 
    xlab("") + 
    ylab("") + 
    theme(axis.text.y = element_text(size=text.size))
  
    print(aa)

    if(save){
       pdf(paste0(name, ".pdf"), height=height, width=width)
       print(aa)
       dev.off()
    }
}


go_barplot_simple <- function(name, dat, save=FALSE, height=12, width=12){
  aa <- dat %>% 
    ggplot(aes(x=Term, y=Score)) + 
    theme_classic() + 
    geom_col(position = "dodge")  + 
    coord_flip() +
    theme(legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_fill_brewer(palette = "Set1") + 
    ggtitle(name) + 
    xlab("") + 
    ylab("-log10 (FDR)") + 
    theme(axis.text.y = element_text(size=text.size))
  
    print(aa)

    if(save){
       pdf(paste0(name, ".pdf"), height=height, width=width)
       print(aa)
       dev.off()
    }
}



### integrate enrichR for multiple groups
# input needs to contains 2 columns
# Gene 
# Group
run_enrichR_groups_multidplyr <- function(gene.list, dbs="GO", prefix="enrichR", core=2,
                                          by.enrich="Combined.Score", top=20, height=10, width=8){
  #require(multidplyr)
  if(dbs == "GO"){
    print(paste("Use preset pathways (GO, KEGG, REACTOME, BIOCARTA)"))
    dbs <- c("GO_Molecular_Function_2017b", "GO_Cellular_Component_2017b", "GO_Biological_Process_2017b" ,
         "BioCarta_2016" ,"KEGG_2016", "Reactome_2016")
  }else if(dbs == "TF"){
    print(paste("Use preset TF pathways"))
    dbs <- c("ChEA_2016", "TRANSFAC_and_JASPAR_PWMs", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
             "ENCODE_TF_ChIP-seq_2015" ,"Transcription_Factor_PPIs")
  }else{
    print(paste("User defined dbs"))
  }
  
  type.col <- gene.list %>% 
    dplyr::select(-Gene) %>% 
    names
  
  enrich.res <- gene.list %>% 
    dplyr::rename(Type = !!sym(type.col)) %>% 
    group_by(Type) %>% 
    nest()
  idx <- rep(1:core, length.out = nrow(enrich.res))
  enrich.res <- bind_cols(tibble(idx), enrich.res)

  by_group <- multidplyr::partition(enrich.res, idx, cluster = create_cluster(core))

  enrich.res <- by_group %>% 
     # Assign libraries
     cluster_library("tidyverse") %>%
     cluster_library("enrichR") %>%
     # Assign values (use this to load functions or data to each core)
     cluster_assign_value("dbs", dbs) %>%
     mutate(res = map(data, ~ enrichr(unlist(.x), dbs))) %>% 
     collect() %>% 
     ungroup() %>% 
     dplyr::select(-idx)
  
  # split by dbs
  by_dbs <- enrich.res %>% 
    dplyr::select(-data) %>% 
    unnest(res) %>% 
    bind_cols(tibble(dbs = rep(dbs, length.out=nrow(.)))) %>%
    split(dbs) %>% 
    map(~ unnest(.x, res) %>% extract_top_enrichR_res(top=top, by="Combined.Score"))
  names(by_dbs) <- substr(names(by_dbs), 1, 30)
  writexl::write_xlsx(by_dbs, paste0(prefix, ".bydbs.xlsx"))
  
  pdf(paste0(prefix, ".top", top, ".pdf"), height=height, width=width)
  for(name in names(by_dbs)){
    cur.dat <- by_dbs[[name]]
    if(nrow(cur.dat) >= 2 & ncol(cur.dat) >= 3){
       x2 <- convert_tibble_to_dataframe(cur.dat)
       fontsize_row = pheatmap_fontsize_row(x2)
       pheatmap(x2, main=name, fontsize_row = fontsize_row, clustering_method = "ward.D2")
    }
  }
  dev.off()

  # split by Type
  by_type <- enrich.res %>% 
    dplyr::select(-data) %>% 
    unnest(res) %>% 
    bind_cols(tibble(dbs = rep(dbs, length.out=nrow(.)))) %>%
    unnest(res) %>% 
    split(.$Type) 
  writexl::write_xlsx(by_type, paste0(prefix, ".byType.xlsx"))

  ### three columns for go_barplot, Term, Score, Pathway
  # min.score <- ifelse(by.enrich == "Adjusted.P.value", 3, 20)
  pdf(paste0(prefix, ".bar.pdf"), height=10, width=12)
  aa <- by_type %>% 
    map(~ .x %>% 
          dplyr::rename(Score = Combined.Score) %>% 
          mutate(Pathway = dbs) %>% 
          dplyr::select(Term, Score, Pathway) %>%
          group_by(Pathway) %>%
          arrange(Pathway, desc(Score)) %>%
          dplyr::slice(1:10) %>%
          arrange(Pathway, Score) %>%
          ungroup() %>%
          dplyr::mutate(Term = gsub(" \\(GO.*", "", Term)) %>%
          dplyr::mutate(Term = gsub("_Homo.*", "", Term)) %>%
          dplyr::mutate(Term = factor(Term, levels=Term %>% unique))
    ) %>% 
    map2(paste(names(.), by.enrich, sep="."), ., go_barplot)
  dev.off()
}




### integrate enrichR for multiple groups
# input needs to contains 2 columns
# Gene 
# Group
run_enrichR_groups <- function(gene.list, dbs="GO", prefix="enrichR", 
                               by.enrich="Combined.Score", top=20, height=10, width=8){
  if(dbs == "GO"){
    print(paste("Use preset pathways (GO, KEGG, REACTOME, BIOCARTA)"))
    dbs <- c("GO_Molecular_Function_2017b", "GO_Cellular_Component_2017b", "GO_Biological_Process_2017b" ,
         "BioCarta_2016" ,"KEGG_2016", "Reactome_2016")
  }else if(dbs == "TF"){
    print(paste("Use preset TF pathways)"))
    dbs <- c("ChEA_2016", "TRANSFAC_and_JASPAR_PWMs", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
             "ENCODE_TF_ChIP-seq_2015" ,"Transcription_Factor_PPIs")
  }else{
    print(paste("User defined dbs"))
  }
  
  type.col <- gene.list %>% 
    dplyr::select(-Gene) %>% 
    names
  enrich.res <- gene.list %>% 
    dplyr::rename(Type = !!sym(type.col)) %>% 
    group_by(Type) %>% 
    nest() %>% 
    mutate(res = map(data, ~ enrichr(unlist(.x), dbs)))
  
  # split by dbs
  by_dbs <- enrich.res %>% 
    dplyr::select(-data) %>% 
    unnest(res) %>% 
    bind_cols(tibble(dbs = rep(dbs, length.out=nrow(.)))) %>%
    split(dbs) %>% 
    map(~ unnest(.x, res) %>% extract_top_enrichR_res(top=top, by="Combined.Score"))
  names(by_dbs) <- substr(names(by_dbs), 1, 30)
  writexl::write_xlsx(by_dbs, paste0(prefix, ".bydbs.xlsx"))
  
  pdf(paste0(prefix, ".top", top, ".pdf"), height=height, width=width)
  for(name in names(by_dbs)){
    cur.dat <- by_dbs[[name]]
    if(nrow(cur.dat) >= 2 & ncol(cur.dat) >= 3){
       x2 <- convert_tibble_to_dataframe(cur.dat)
       fontsize_row = pheatmap_fontsize_row(x2)
       pheatmap(x2, main=name, fontsize_row = fontsize_row, clustering_method = "ward.D2")
    }
  }
  dev.off()

  # split by Type
  by_type <- enrich.res %>% 
    dplyr::select(-data) %>% 
    unnest(res) %>% 
    bind_cols(tibble(dbs = rep(dbs, length.out=nrow(.)))) %>%
    unnest(res) %>% 
    split(.$Type) 
  writexl::write_xlsx(by_type, paste0(prefix, ".byType.xlsx"))

  ### three columns for go_barplot, Term, Score, Pathway
  # min.score <- ifelse(by.enrich == "Adjusted.P.value", 3, 20)
  pdf(paste0(prefix, ".bar.pdf"), height=10, width=12)
  aa <- by_type %>% 
    map(~ .x %>% 
          dplyr::rename(Score = Combined.Score) %>% 
          mutate(Pathway = dbs) %>% 
          dplyr::select(Term, Score, Pathway) %>%
          group_by(Pathway) %>%
          arrange(Pathway, desc(Score)) %>%
          dplyr::slice(1:10) %>%
          arrange(Pathway, Score) %>%
          ungroup() %>%
          dplyr::mutate(Term = gsub(" \\(GO.*", "", Term)) %>%
          dplyr::mutate(Term = gsub("_Homo.*", "", Term)) %>%
          dplyr::mutate(Term = factor(Term, levels=Term %>% unique))
    ) %>% 
    map2(paste(names(.), by.enrich, sep="."), ., go_barplot)
  dev.off()
}






# parse/summarise count data 
summarise_cntdata <- function(dat, gene_col="Gene_name", 
                              left_col="ID", right_col="chromosome_name"){
   left.idx <- which(colnames(dat) == left_col)
   right.idx <- which(colnames(dat) == right_col)
   
   res <- dat %>% 
     dplyr::select(-(1:left.idx), -(right.idx:ncol(dat))) %>% 
     dplyr::select(!!rlang::sym(gene_col), everything()) %>% 
     group_by(!!rlang::sym(gene_col)) %>% 
     summarise_all(sum)
}

check_ctrl_sgRNA <- function(rawdata, pattern="_IGO", ctrl="NonTargeting", prefix="ctrl.corr"){
   require(GGally)

   controlData <- rawdata %>% 
      dplyr::select(-sgRNA_Target_Sequence) %>% 
      set_names(gsub(paste0(pattern, ".*"), "", colnames(.))) %>% 
      filter(rowSums(.[, -1]) > 0) %>% 
      dplyr::mutate(sgRNA = sgRNA_ID) %>% 
      mutate(sgRNA = gsub("[-_].*", "", sgRNA_ID)) %>% 
      mutate(sgRNA = ifelse(grepl(ctrl, sgRNA), "Control", "sgRNA"))

   pdf(paste0(prefix, ".normal.pdf"), height=12, width=12)
   print(GGally::ggpairs(controlData, 
                   columns = 2:ncol(controlData),
                   # columns = 2:4,
                   mapping = ggplot2::aes(colour=sgRNA),
                   # lower=list(continuous='smooth')),
                   upper = list(continuous = function(data, mapping, ...) {
                                   ggally_cor(data = data, mapping = mapping, alignPercent = 0.9, size=4) + 
                                      scale_colour_brewer(palette = "Set1") + 
                                      theme_classic()}),
                   lower = list(continuous = function(data, mapping, ...) {
                                   ggally_smooth(data = data, mapping = mapping, alpha = 0.6) + 
                                      scale_colour_brewer(palette = "Set1") + 
                                      theme_bw()}),
                   diag = list(continuous = function(data, mapping, ...) {
                                  ggally_densityDiag(data = data, mapping = mapping, alpha = 0.7) + 
                                     scale_fill_brewer(palette = "Set1") + 
                                     theme_bw()})
                   ) +  theme(strip.text = element_text(size=14)))
   dev.off()
   
   pdf(paste0(prefix, ".log.pdf"), height=12, width=12)
   print(GGally::ggpairs(controlData, 
                   columns = 2:ncol(controlData),
                   # columns = 2:4,
                   mapping = ggplot2::aes(colour=sgRNA),
                   # lower=list(continuous='smooth')),
                   upper = list(continuous = function(data, mapping, ...) {
                                   ggally_cor(data = data, mapping = mapping, alignPercent = 0.9, size=4) + 
                                      scale_colour_brewer(palette = "Set1") + 
                                      theme_classic()}),
                   lower = list(continuous = function(data, mapping, ...) {
                                   ggally_smooth(data = data, mapping = mapping, alpha = 0.6) + 
                                      scale_x_log10() + 
                                      scale_y_log10() + 
                                      scale_colour_brewer(palette = "Set1") + 
                                      theme_bw()}),
                   diag = list(continuous = function(data, mapping, ...) {
                                  ggally_densityDiag(data = data, mapping = mapping, alpha = 0.7) + 
                                     scale_fill_brewer(palette = "Set1") + 
                                     scale_x_log10() + 
                                     scale_y_log10() + 
                                     theme_bw()})
                   ) +  theme(strip.text = element_text(size=14)))
   dev.off()
}


# for genomeCRISPR mutation_coloring 
merge_mutation <- function(mutation, type = "simple"){
   if(tolower(type) == "simple"){
      mutation <- case_when(
                            mutation == "WT" ~ "WT",
                            mutation == "Silent" ~ "WT",
                            mutation == "NODATA" ~ "NOINFO",
                            mutation == "NOCCLE" ~ "NOINFO",
                            TRUE ~ "Altered"
                            )
   }else if(tolower(type) == "merge"){
      mutation <- case_when(
                            mutation == "WT" ~ "WT",
                            grepl("Silent", mutation) ~ "WT",
                            grepl("In_Frame", mutation) ~ "WT",
                            mutation == "NODATA" ~ "NOINFO",
                            mutation == "NOCCLE" ~ "NOINFO",
                            grepl("Missense", mutation) ~ "Missense",
                            grepl("Nonsense", mutation) ~ "Nonsense",
                            grepl("Splice_Site", mutation) ~ "Splice_Site",
                            grepl("Frame_Shift", mutation) ~ "Splice_Site",
                            TRUE ~ "Altered"
                            )
   }else if(tolower(type) == "all"){
      mutation <- case_when(
                            mutation == "WT" ~ "WT",
                            grepl("Silent", mutation) ~ "Silent",
                            grepl("In_Frame", mutation) ~ "In_Frame",
                            mutation == "NODATA" ~ "NODATA",
                            mutation == "NONCCLE" ~ "NONCCLE",
                            grepl("Missense", mutation) ~ "Missense",
                            grepl("Nonsense", mutation) ~ "Nonsense",
                            grepl("Splice_Site", mutation) ~ "Splice_Site",
                            grepl("Frame_Shift", mutation) ~ "Frame_Shift",
                            TRUE ~ "Others"
                            )
   }else{
      print(paste("### this type:", type, "cannot be recognized, try again")) 
      return(NA)
   }
   return(mutation)
}

mut_color_predefined <- function(name = "c1"){
   if(name == "c1"){
      my_pallete <- c(
                      "WT" = rgb(204,204,204, maxColorValue=255),
                      "Missense" = rgb(56,146,208, maxColorValue=255),
                      "Altered" = rgb(255,140,10, maxColorValue=255),
                      "Splice_Site" = rgb(116,187,88, maxColorValue=255),
                      "Others" = rgb(236,67,35, maxColorValue=255),
                      "Silent" = rgb(205,95,46, maxColorValue=255),
                      "In_Frame" = rgb(125,95,46, maxColorValue=255),
                      "Nonsense" = rgb(100,100,100, maxColorValue=255),
                      "NOINFO" = rgb(203,77,202, maxColorValue=255),
                      "NODATA" = rgb(203,77,202, maxColorValue=255),
                      "NONCCLE" = rgb(203,77,202, maxColorValue=255),
                      "Frame_Shift" = rgb(115,113,206, maxColorValue=255))
      return(my_pallete)
   }else{
      print(paste("### Currently only support name: `c1`, please try again"))
      return(NA)
   }
}

pdac_color_predefined <- function(name = "paper"){
   if(name == "paper"){
                      #more.samples <- c("T1.Normal", "T2.Reg-ADM", "T3.KrasG12D", "T4.Tum-ADR", "T5.PanIN-early", "T6.PanIN-late", "T7.PDAC")
                      #more.color <- c("gray34", "green4", "brown", "orange2", "darkorange4", "mediumpurple2", "blue4")
      my_pallete <- c(
                      #"T1.Normal" = rgb(87,87,87, maxColorValue=255),
                      "T1.Normal" = rgb(147,147,147, maxColorValue=255),
                      #"T2.Reg-ADM" = rgb(0,139,0, maxColorValue=255),
                      "T2.Reg-ADM" = rgb(0,169,0, maxColorValue=255),
                      #"T3.KrasG12D" = rgb(165,42,42, maxColorValue=255), 
                      "T3.KrasG12D" = rgb(225,42,42, maxColorValue=255), 
                      "T4.Tum-ADR" = rgb(238,154,0, maxColorValue=255),
                      "T5.PanIN-early" = rgb(139,69,0, maxColorValue=255),
                      "T6.PanIN-late" = rgb(159,121,238, maxColorValue=255),
                      "T7.PDAC" = rgb(0,0,139, maxColorValue=255),
                      "T8.PDAC-2D" = rgb(130,20,139, maxColorValue=255), 
                      "T10.Normal-shRen" = rgb(87,87,87, maxColorValue=255), 
                      "T11.Normal-shRen" = rgb(100,200,32, maxColorValue=255), 
                      "T12.Reg-ADM-shRen" = rgb(0,139,0, maxColorValue=255), 
                      "T13.Reg-ADM-sh552" = rgb(130,20,39, maxColorValue=255), 
                      "T14.Reg-ADM-sh1448" = rgb(128,255,204, maxColorValue=255), 
                      "T15.Tum-ADR-shRen" = rgb(238,154,0, maxColorValue=255),
                      "T16.Tum-ADR-sh552" = rgb(100,210,109, maxColorValue=255), 
                      "T17.Tum-ADR-sh1448" = rgb(255,212,128, maxColorValue=255)) 
      return(my_pallete)
   }else{
      print(paste("### Currently only support name: `paper`, please try again"))
      return(NA)
   }
}

### from plotflow
# https://github.com/trinker/plotflow/blob/master/R/mergePDF.R
mergePDF <-
  function(..., file, gsversion = NULL, in.file = NULL) {
    if (is.null(in.file)) {
      in.file <- substitute(...())
    }
    infiles <- unlist(lapply(in.file, function(y) as.character(y)))
    if (is.null(gsversion)) {
      gsversion <- tools::find_gs_cmd()
      if (!nzchar(gsversion))
        stop("Please install Ghostscript and see ?tools::find_gs_cmd")
    }
    pre <- c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite")
    out <- paste0("-sOutputFile=", shQuote(file))
    # system2(gsversion, c(pre, out, shQuote(infiles)))
    system2(gsversion, c(pre, out, infiles))
}


### boxplot
pdac_boxplot <- function(dat, Gene_name="GeneA", prefix = "all", height=10, width=20){
  p1 <- dat %>%
    ggplot(aes(x=Stage_allinfo, y=value, color=Stage_allinfo)) +
    geom_boxplot(outlier.colour = NULL, outlier.alpha = 0) +
    geom_jitter(width=0.2, alpha=0.6) +
    theme_bw() +
    xlab("") +
    ylab("") +
    ggtitle(Gene_name) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=24, face="bold"),
          axis.text.x = element_text(angle=90, vjust=0.5, size=14, hjust=0.95), 
          axis.text.y = element_text(size=14), 
          legend.key.size = unit(2, 'lines'),
          legend.key = element_rect(size = 8, color = NA),
          legend.text = element_text(size = 12),
          legend.title = element_blank()) + 
    guides(fill = guide_legend(nrow=4, byrow=TRUE), color = guide_legend(nrow=4, byrow=TRUE))  
  
  pdf(paste(prefix, Gene_name, 'pdf', sep="."), height=height, width=width)
    print(p1)
  dev.off()
}

pdac_boxplot_prog <- function(dat.select, Gene_name="GeneA", color = "black", prefix = "t1-4", height=10, width=8){
   p2 <- dat.select %>%
    ggplot(aes(x=Stage_allinfo, y=value, fill=Stage_allinfo, color=Stage_allinfo)) +
    geom_boxplot(outlier.colour = NULL, outlier.alpha = 0, alpha=0.4) +
    geom_jitter(width=0.2, alpha=0.6) +
    theme_bw() +
    xlab("") +
    ylab("") +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    ggtitle(Gene_name) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=24, face="bold"),
          axis.text.x = element_text(angle=90, vjust=0.5, size=14, hjust=0.95), 
          axis.text.y = element_text(size=14), 
          legend.key.size = unit(2, 'lines'),
          legend.key = element_rect(size = 8, color = NA),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) + 
    guides(fill = guide_legend(nrow=2, byrow=TRUE), color = guide_legend(nrow=2, byrow=TRUE))  

  pdf(paste(prefix, Gene_name, 'pdf', sep="."), height=height, width=width)
    print(p2)
  dev.off()
}


# Extract the feature Genes in each group
wgcna_extract_hubgenes <- function(datExpr, MM.cutoff=0.7, math.cutoff.quantile=3, math="mad", top.math=1500, plot=F, save.data=F, prefix="MM_vs_MAD"){
   require(WGCNA)
   require(ggplot2)

   MM.cutoff <- MM.cutoff
   math.Expr <- log(apply(datExpr, 2, math))
   Exp.cutoff <- signif(quantile(math.Expr)[math.cutoff.quantile], digits = 3)

   mad.genes <- extract_data_by_mad(datExpr, topN=top.math, by="col", type="genes")

   # OVData.top1500mad.genes <- rownames(top.mad.data)

   total.num <- 0
   total.genes <- data.frame(Gene=NULL, Cor=NULL, Module=NULL)

   for(module in modNames){
     column = match(module, modNames);
     moduleGenes = moduleColors==module;
     current.module <- data.frame(Gene=colnames(datExpr)[moduleGenes],
                                   Exp=math.Expr[moduleGenes],
                                   MM=abs(geneModuleMembership[moduleGenes, column]),
                                   Module=module)

     mad.module <- subset(current.module, Gene %in% mad.genes)
     Cor <- signif(cor(current.module$Exp, current.module$MM), 2)
     Corp <- signif(corPvalueStudent(Cor, sum(is.finite(current.module$Exp) & is.finite(current.module$MM))), 2)

     if(plot){
         a <- ggplot(current.module, aes(x=MM, y=Exp, label=Gene)) +
               geom_point(alpha=0.5, colour=module) +
               geom_text(size=1.5, hjust = 0, nudge_x = 0.005) + theme_bw() +
               geom_vline(xintercept = MM.cutoff, colour="tomato") +
               geom_hline(yintercept = Exp.cutoff, colour="tomato") +
               xlab(paste("Module Membership in", module, "module")) +
               ylab("Log (expression)") +
               ggtitle(paste("Log expression vs. Module membership\n", module, "module:", sum(moduleGenes), "genes\n", "cor=", Cor, "p=", Corp)) +
               theme(plot.title = element_text(size=16, face="bold")) +
               theme(axis.text=element_text(size=16), axis.title=element_text(size=16, face="bold")) +
               theme(strip.text.x = element_text(size=16, face="bold")) +
               theme(axis.title.y=element_text(margin=margin(0,20,0,10))) +
               theme(axis.title.x=element_text(margin=margin(20,0,10,0)))

         # add red dots only if we have genes that are within TopMAD
         if(nrow(mad.module)>0){
           a <- a + geom_point(alpha=0.7, aes(x=MM, y=Exp), colour="red", data=mad.module) +
                  geom_text(size=1.5, hjust = 0, nudge_x = 0.005, data=mad.module, colour="red") 
         }
         (gg <- ggplotly(a))
         # print(gg)
         htmltools::tagList(gg)
      }

      idx <- (abs(current.module$MM) >= MM.cutoff) & (current.module$Exp >= Exp.cutoff)
      total.num <- total.num + sum(idx)
      total.genes <- rbind(total.genes, current.module[idx,])
   }

   print(paste("total", total.num, "genes have cor above", MM.cutoff, "and Log(Expression) above", Exp.cutoff))
   if(save.data){
      # filename <- paste0(cancer, "_TopMAD", topN, "_Module_Genes_MMcutoff", MM.cutoff, "_Expcutoff", Exp.cutoff, ".txt")
      filename <- paste0(prefix, "_Module_Genes_MMcutoff", MM.cutoff, "_Expcutoff", Exp.cutoff, ".txt")
      write.table(total.genes, filename, sep="\t", quote=F, row.names=F)
   }

   return(total.genes)
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



assign_dir <- function(x){
  res <- x %>% 
    dplyr::select(Gene = Gene_name, log2FC = log2FoldChange) %>% 
    group_by(Gene) %>% 
    filter(abs(log2FC) == max(abs(log2FC))) %>% 
    ungroup(Gene) %>% 
    mutate(Group = ifelse(log2FC > 0, "UP", "DN")) 
  return(res)
}

# input list of two DESeq2 results 
# list(dat1 = xxx, dat2 = yyy)
# dat1 and dat2 will be used as filename of the final table
# assume it's house, if convert = TRUE, convert to Human
compare_de_res <- function(x, convert=FALSE){
  filename <- names(x)
  comp_res <- tibble(Sample = filename, content = x) %>%
    mutate(DE = map(content, assign_dir)) %>% 
    dplyr::select(Sample, DE) %>% 
    unnest(DE) %>% 
    unite(Sample, Sample, Group) 
 
  if(convert){
     comp_res <- comp_res %>% 
        convert_mouse_to_human() 
  }
  comp_res <- comp_res %>% 
        group_by(Gene, Sample) %>% 
        filter(abs(log2FC) == max(abs(log2FC))) %>% 
        dplyr::select(-log2FC) %>% 
        ungroup() %>% 
        mutate(Dir = gsub(".*_", "", Sample)) 

  up <- comp_res %>%  
    filter(Dir == "UP") %>% 
    group_by(Gene) %>% 
    summarise(cnt = n()) %>% 
    left_join(., comp_res)  %>% 
    filter(!is.na(Gene), Dir == "UP") %>% 
    mutate(Sample = ifelse(cnt == 2, paste0(paste0(filename, collapse = "_"), "_UP_shared"), Sample)) %>% 
    distinct(Gene, Sample) 

  dn <- comp_res %>%  
    filter(Dir == "DN") %>% 
    group_by(Gene) %>% 
    summarise(cnt = n()) %>% 
    left_join(., comp_res)  %>% 
    filter(!is.na(Gene), Dir == "DN") %>% 
    mutate(Sample = ifelse(cnt == 2, paste0(paste0(filename, collapse = "_"), "_DN_shared"), Sample)) %>% 
    distinct(Gene, Sample) 
  
  final_res <- bind_rows(up, dn) %>% 
    group_by(Sample) %>% 
    mutate(idx = 1:n()) %>% 
    ungroup() %>% 
    mutate(Sample = factor(Sample, levels=unique(Sample))) %>% 
    spread(Sample, Gene) %>% 
    replace(is.na(.), "") %>% 
    dplyr::select(-idx)
  
  return(final_res)
}


my_DESeq2_pcaplot <- function(object, intgroup = "condition", genes = NULL, returnData = FALSE){

   #select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
   pca <- prcomp(t(object[rownames(object) %in% genes, ]))
   percentVar <- pca$sdev^2/sum(pca$sdev^2)

   if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
   }
   intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                drop = FALSE])
   group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
   }else {
      colData(object)[[intgroup]]
   }
   d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                   intgroup.df, name = colnames(object))
   if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
   }
   ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
      coord_fixed()
}


### 
compare_group <- function (group1, group2, file.prefix = "Compare.group", title = "", 
    plot = T, save.image = F, add.legend = T, label.size = 5, 
    title.size = 36, angle = 45, save.data=F) 
{
    if (class(group1) == "integer" | class(group1) == "numeric") 
        group1 <- paste0("Group", group1)
    if (class(group2) == "integer" | class(group2) == "numeric") 
        group2 <- paste0("Group", group2)
    group1 <- factor(group1)
    group2 <- factor(group2, levels=paste0("A", c(4,3,1,2,5,6, "O")))
    summary <- as.matrix(table(group1, group2))
    sums <- rowSums(summary)
    summary.pct <- apply(t(summary), 1, function(x) x/sums) * 
        100
    summary.pct <- melt(summary.pct)
    if(save.data){
       write_csv(summary.pct, paste0(file.prefix, ".csv"))
    }

    colourCount = length(unique(summary.pct$group2))
    num <- ifelse(colourCount > 9, 9, colourCount)
    getPalette = colorRampPalette(RColorBrewer::brewer.pal(num, 
        "Set1"))
    ray <- ggplot(data = summary.pct, aes(x = group1, y = value, 
        fill = group2)) + geom_bar(stat = "identity") + labs(y = "Percentage (%)") + 
        #scale_fill_manual(values = getPalette(colourCount)) + 
        #scale_fill_manual(values = c("#F4C6C4", "#B8F9D6", "#FF836F", "#941100", "#76D6FF", "#011993", "gray80", "green")) + 
        scale_fill_manual(values = c("#941100", "#FF836F", "#F4C6C4", "#B8F9D6", "#76D6FF", "#011993", "gray80", "green")) + 
        xlab("") + labs(fill = "")
    naikai <- ggplot_build(ray)$data[[1]]
    naikai$position = (naikai$ymax + naikai$ymin)/2
    foo <- summary.pct %>% arrange(group1, desc(group2)) %>% 
        .$value %>% round(digits = 2)
    naikai$label <- paste0(foo, "%")
    ray <- ray + annotate(x = naikai$x, y = naikai$position, 
        label = naikai$label, geom = "text", size = label.size)
    ray <- ray + ggtitle(title) + theme(plot.title = element_text(face = "bold", 
        size = title.size))
    ray <- ray + theme(axis.text.x = element_text(angle = angle, 
        vjust = 1, hjust = 1))
    if (!add.legend) {
        ray <- ray + theme(legend.position = "none")
    }
    if (plot) 
        print(ray)
    if (save.image) {
        pdf(paste0(file.prefix, ".pdf"), height = 10, width = 12)
        print(ray)
        dev.off()
    }
    return(ray)
}

blank_theme <- theme_minimal() +
   theme(
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.border = element_blank(),
         panel.grid=element_blank(),
         axis.ticks = element_blank(),
         plot.title=element_text(size=14, face="bold")
         )

deseq2gseaRnk <- function(deseq, gene="Gene", log2fc="log2FoldChange", pvalue="pvalue", 
                          min.exp=0, species="mouse", type="signp"){
   if(tolower(type) == "signp"){
      deseq.rnk <- deseq %>% 
         mutate(Rank = sign(!!sym(log2fc)) * -log10(!!sym(pvalue))) 
   }else if(tolower(type) == "log2fc"){
      deseq.rnk <- deseq %>% 
         mutate(Rank = !!sym(log2fc)) 
   }else{
      warnings(paste("Unknown type:", type))
   }

  deseq.rnk <- deseq.rnk %>% 
    filter(!is.na(Rank)) %>%
    group_by(!!sym(gene_col)) %>% 
    arrange(desc(Rank)) %>% 
    dplyr::slice(1) %>% 
    ungroup() %>% 
    arrange(desc(Rank)) %>% 
    filter(baseMean >= min.exp) 

   if(toupper(species) == "MOUSE"){
      rnk <- deseq.rnk %>% 
       dplyr::select(Gene = !!sym(gene), Rank) %>%
       convert_human_to_mouse()
   }else if(toupper(species) == "HUMAN"){
      rnk <- deseq.rnk %>% 
        dplyr::select(Gene = !!sym(gene), Rank) %>%
        convert_mouse_to_human()
   }else{
      warnings(paste("Unknown species:", species))
   }

   return(rnk)
}
