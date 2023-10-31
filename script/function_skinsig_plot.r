

## Required input:
# - count_mx = count matrix (gene symbol in row.names and duplicates are not allowed; missing values not allowed in the count matrix, unique sample id in the column)
# - group = vector indicating test or control for the samples (must be in the same order as the count matrix); can only do 2 group comparison
# - ctrl_name = the value in group that indicates control samples (default = "ctrl")
# - datatype = the data type (default = "RNAseq")
# - proportion_same_direction = the proportion of genes regulated in the same direction in gene set enrichment analysis (default = 0.75)
# - proportion_same_direction = FDR threshold for gene set enrichment analysis (default = 0.01)

skinsig_plots <- function(count_mx, group, ctrl_name = "ctrl", datatype="RNAseq", proportion_same_direction = 0.75, fdr_threshold = 0.01){
	## declare outputs
	out <- list()
	test_name <- unique(group[group != ctrl_name])
	group <- factor(group, levels=c(ctrl_name, test_name))
	
	## Data processing
	# make DEGlist
	deglist <- edgeR::DGEList(counts = count_mx, group = group)
	# remove lowly expressed genes
	keep.exprs <- edgeR::filterByExpr(deglist, group=group)
	#keep.exprs[names(keep.exprs) %in% skinsig$gene] <- TRUE # keep all genes that are present in skinsig
	deglist <- deglist[keep.exprs,, keep.lib.sizes=FALSE]
	# normalise
	deglist <- edgeR::calcNormFactors(deglist, method = "TMM")
	#deglist <- DESeq2::estimateDispersionsFit(deglist)
	# fit stats over group
	design <- model.matrix(~group)
	if(datatype == "RNAseq"){
		deglist2 <- voom(deglist, design)
	} else {
		deglist2 <- deglist
	}
	vfit <- limma::lmFit(deglist2, design)
	efit <- limma::eBayes(vfit)
	top_genes <- limma::topTable(efit, n=Inf, sort.by="none")
	
	## Gene set enrichment with mroast
	# Use mroast to determine if each skinsig is significant
	# Identify corresponding index for each signature
	current_skinsig <- skinsig[skinsig$hgnc_symbol %in% row.names(deglist),]
	skinsig_list <- split(current_skinsig, current_skinsig$signature)
	idx_match <- lapply(skinsig_list, function(x)row.names(deglist) %in% x$hgnc_symbol)
	set.seed(123) # random process in mroast
	mroast_result <- limma::mroast(deglist2, index=idx_match, design, contrast=ncol(design))
	out[["enrichment"]] <- mroast_result
	
	## Plot heatmap
	# Find the logFC for all genes
	top_genes <- top_genes[row.names(deglist2),] # make sure the order of genes is in the original order
	logfc <- lapply(idx_match, function(x)mean(top_genes$logFC[x]))
	logfc <- do.call(c, logfc)
	plot_df <- data.frame(logFC = logfc, skinsig = names(logfc))
	plot_df <- merge(plot_df, mroast_result, by.x="skinsig", by.y=0)
	plot_df$skinsig_significant <- ifelse((plot_df$PropDown >= proportion_same_direction | plot_df$PropUp >= proportion_same_direction ) & (plot_df$FDR <= fdr_threshold), 1, 0)
	# ggplot geom_tile()
	p <- ggplot(data= plot_df, aes(x=skinsig, y = 1, fill=logFC)) + 
		geom_tile() + 
		geom_point(aes(alpha=skinsig_significant), colour="#ffffff") +
		coord_fixed(ratio = 1) + 
		theme_bw() + xlab("") + ylab("") + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y =element_blank(), axis.ticks.y=element_blank()) +
		scale_alpha_continuous(range = c(0, 1)) +
		scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint =0) +
		guides(alpha = "none")

	out[["heatmap"]] <- p
	
	
	## Plot within-signature correlation
	lcmp <- edgeR::cpm(deglist2)
	correlation_skinsig <- lapply(idx_match, function(x){corr_mx <- cor(lcmp[x,]); corr_mx[lower.tri(corr_mx, diag = TRUE)] <- NA; return(as.vector(corr_mx))})
	plot_df_corr <- lapply(names(correlation_skinsig), function(x)data.frame(skinsig = x, correlation = correlation_skinsig[[x]]))
	plot_df_corr <- do.call(rbind, plot_df_corr)
	plot_df_corr <- plot_df_corr[!is.na(plot_df_corr$correlation), ]
	p2 <- ggplot(data=plot_df_corr, aes(x=skinsig, y=correlation)) + geom_jitter(size=0.1) + 
		theme_bw() + xlab("") +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	print(p2)
	print(p)

	out[["correlation"]] <- p2

	return(out)

	## 
}