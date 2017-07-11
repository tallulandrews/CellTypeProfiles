# Generalizable - Step 1 = Get profiles
get_cluster_profile <- function(expr_mat, clusters, norm_method="CPM", is.log=FALSE, out.log=2, feature_selection=kw.features) {
	if (min(factor_counts(clusters)) < 2) {stop("Error: Cannot have a cell-type with a single sample.")}
	# is.log should hold the base of the log if it has been log-transformed.
	# If no norm here, need at least FS!
	if (is.log) {
		expr_mat <- is.log^expr_mat;
		expr_mat <- expr_mat - 1;
	}
	profiles <- my_row_mean_aggregate(expr_mat, clusters)
	is.feature <- feature_selection(expr_mat, clusters)
	if (norm_method=="CPM") {
		# CPM
		normfact <- colSums(profiles);
	} else if (norm_method %in% c("TMM", "RLE", "upperquartile", "none")) {
		# TMM
		#require("edgeR")
		normfact <- edgeR::calcNormFactors(as.matrix(profiles), method=norm_method)
	} else {
		warning("Warning: Unrecognized normalization method. Profiles will be scaled but not normalized.")
		normfact <- rep(1, times=length(profiles[1,]))
	}

	# Apply & scale to ~ 1,000,000 counts per profile
	profiles <- t(t(profiles)/normfact)
	scalefact <- 1000000/median(colSums(profiles))
	profiles <- profiles*scalefact

	 # out.log should hold the base of the log if should be log-transformed.
	if (out.log) {
		profiles <- log(profiles+1)/log(out.log);
	}

	return(cbind(profiles, is.feature))
}

# Generalizable - Step 2 = Match clusters using reciprocal best-hits
get_reciprocal_hits <- function(profiles1, profiles2, features=NA) {
	# features = gene names of features to use
        
        dim1 = dim(profiles1);
        dim2 = dim(profiles2);


        if (dim1[1] != dim2[1]) {
                profiles1 <- profiles1[rownames(profiles1) %in% rownames(profiles2),]
                profiles2 <- profiles2[match(rownames(profiles1), rownames(profiles2)),]
                dim1 = dim(profiles1);
                dim2 = dim(profiles2);
        }
        features1 = profiles1[,dim1[2]]
        features2 = profiles2[,dim2[2]]
        profiles1 <- as.matrix(profiles1[,-dim1[2]])
        profiles2 <- as.matrix(profiles2[,-dim2[2]])

        if (is.na(features)) { # throw warning when features is not NA
                # feature of one or other dataset.
                consensus_features = features1==1 | features2==1;
        } else {
                # user provided features.
                consensus_features = rownames(profiles1) %in% features;
        }

        # Get Matches
        match_cors <- cor(profiles1[consensus_features,], profiles2[consensus_features,], method="spearman");
        matches1 = apply(match_cors, 1, function(x){which(x==max(x))})
        matches2 = apply(match_cors, 2, function(x){which(x==max(x))})

        # Get Reciprocal Best-Hits
        recip = which(1:length(matches1) == matches2[matches1])
        recip_table = cbind(recip,matches1[recip])
        colnames(recip_table) = c("profiles1", "profiles2")
#        match_table = cbind(matches1, matches2) # not necessarily of equal length
#        colnames(match_table) = c("profiles1->profiles2", "profiles2->profiles1")
        return(list(matches1=matches1, matches2=matches2, recip=recip_table))
}

combine_and_match_clusters <- function(profile_list, features=NA){
	# features = gene names of features to use
	if (length(profile_list) < 2) {stop("Error: Cannot combine & match profiles from a single dataset")}

	# Dataset names
	if (is.null(names(profile_list))) {
		warning("Warning: no dataset names, creating default names");
		names(profile_list) <- paste("Data", 1:length(profile_list));
	}
	# attach dataset names to cluster names
	for (i in 1:length(profile_list)) {
		colnames(profile_list[[i]]) <- paste(names(profile_list)[i], colnames(profile_list[[i]]), sep="_");
		profile_list[[i]] <-  profile_list[[i]][order(rownames(profile_list[[i]])),]
	}

	
	# for each pair of datasets get reciprocal best-hits
	pairs = combn(names(profile_list), 2)
	RECIP <- vector()
	for(j in 1:length(pairs[1,])) {
		out = get_reciprocal_hits(profile_list[[pairs[1,j]]], profile_list[[pairs[2,j]]], features=features)
		recip = as.data.frame(out$recip)
		recip[,1] <- colnames(profile_list[[pairs[1,j]]])[recip[,1]]
		recip[,2] <- colnames(profile_list[[pairs[2,j]]])[recip[,2]]
		RECIP <- rbind(RECIP, recip)
	}
	# match them all using graph connected components
	#require("igraph")
	matched_clusters <- igraph::components(igraph::graph_from_edgelist(as.matrix(RECIP), directed=FALSE))$membership
	
	# Combine profiles together
	# consensus genes
	consensus_genes = rownames(profile_list[[1]])
	for (i in 2:length(profile_list)) {
		consensus_genes <- consensus_genes[consensus_genes %in% rownames(profile_list[[i]])]	
	}
	# combine profiles
	# duplicate code
	these_profiles <- as.matrix(profile_list[[1]][rownames(profile_list[[1]]) %in% consensus_genes,-dim(profile_list[[1]])[2]])
	colnames(these_profiles) <- colnames(profile_list[[1]])[-dim(profile_list[[1]])[2]]
	Combined <- these_profiles
	Combined_features <- profile_list[[1]][rownames(profile_list[[1]]) %in% consensus_genes,dim(profile_list[[1]])[2]]
	# end duplicate code
	for (i in 2:length(profile_list)) {
		# duplicate code
		these_profiles <- as.matrix(profile_list[[i]][rownames(profile_list[[i]]) %in% consensus_genes,-dim(profile_list[[i]])[2]])
		colnames(these_profiles) <- colnames(profile_list[[i]])[-dim(profile_list[[i]])[2]]
		Combined <- cbind(Combined, these_profiles)
		Combined_features <- cbind(Combined_features, profile_list[[i]][rownames(profile_list[[i]]) %in% consensus_genes,dim(profile_list[[i]])[2]])
		#end duplicate code
	}
	colnames(Combined_features) <- names(profile_list)

	# merge match clusters
	Combined_Clusters <- colnames(Combined)
	for (comp in 1:max(matched_clusters)) {
		Combined_Clusters[Combined_Clusters %in% names(matched_clusters)[matched_clusters==comp]] <- paste("Match",comp, sep="_")
	}
	# dataset labels
	Combined_Dataset <- vector();
	for (i in 1:length(profile_list)) {
		labs <- rep(names(profile_list)[i], times=dim(profile_list[[i]])[2]-1)
		Combined_Dataset <- c(Combined_Dataset, labs)
	}

	return(list(profiles = Combined, labels = Combined_Clusters, dataset = Combined_Dataset, features = Combined_features))
}

# Generalizable - Step 3 = GLM to calculate batch effects 
#		& Step 4 = remove batch effects from all profiles

glm_of_matches <- function(matches) { #This is not optimized
	matched <- factor_counts(matches$labels);
	matched <- names(matched[matched >1])

	matched_matrix <- matches$profiles[,matches$labels %in% matched]
	matched_type <- matches$labels[matches$labels %in% matched]
	matched_batch <- matches$dataset[matches$labels %in% matched]

	batch_lvl = levels(factor(matches$dataset))
	type_lvl = levels(factor(matched_type))
	all_matrix <- matches$profiles
	if (length(unique(matched_batch)) < length(unique(matches$dataset))
		& length(factor_counts(matches$labels)) > 1) {
		# Insufficient matched groups
		# and Not all groups matched to a single group
		effect_vec = rep(1, times=length(matched_batch));
		all_effect <- rep(1, times=length(matches$dataset))

		for (i in 2:length(batch_lvl)) {
			effect_vec[matched_batch == batch_lvl[i]] <- i;
			all_effect[matches$dataset == batch_lvl[i]] <- i;
		}	
		warning("Warning: Matching groups do not span all datasets, group ID will not be included in GLM.")
		calc_batch_effect_naive <- function(x) {
			effects = c(0,glm(x~matches$dataset)$coef); # Do I force the intercept to be zero i.e. ~0+matches$dataset
			return(effects)
		}
		batch_effects <- apply( matches$profiles, 1, calc_batch_effect_naive)
		corrected <- all_matrix - t(batch_effects[all_effect,])
	} else {
		# Matched groups cover all batches
		effect_vec = rep(1, times=length(matched_batch));
		all_effect <- rep(1, times=length(matches$dataset))
	
		for (i in 2:length(batch_lvl)) {
			effect_vec[matched_batch == batch_lvl[i]] <- i+length(type_lvl);
			all_effect[matches$dataset == batch_lvl[i]] <- i+length(type_lvl);
		}	
		calc_batch_effect <- function(x) {
			effects = c(0,glm(x~matched_type+matched_batch)$coef); # Do I force the intercept to be zero i.e. ~0+matched_type+matched_batch
			return(effects)
		}
		batch_effects <- apply(matched_matrix, 1, calc_batch_effect)
		corrected <- all_matrix - t(batch_effects[all_effect,])
	}
	return(list(corrected_mat=corrected, model_effects=batch_effects))
}

# Generalizable - Step 5 = Cluster corrected profiles
#		& Step 6 = Calculate meta-profiles
cluster_profile_heatmap <- function(corrected_mat, matches, features_only=TRUE, show_genes=FALSE, npermute=0, distfun=function(x){as.dist(1-cor(t(x), method="spearman"))}, hclustfun=function(x){hclust(x,method="complete")}, ann=NULL) {
	#source("/nfs/users/nfs_t/ta6/R-Scripts/heatmap.3.R")

	my_profiles <- corrected_mat
	# Constrain to union of features (optional)
	if (features_only) {
		is.feature = rowSums(matches$features) > 0
		my_profiles <- my_profiles[is.feature,]
	}

	# Permutation significance
	D <- vector();
	M <- my_profiles;
	ncol <- length(my_profiles[1,])
	nrow <- length(my_profiles[,1])
	pairwise_threshold = 0
	if (npermute > 0) {
		set.seed(101); # reproduciblity
		for(rep in 1:npermute) {
			P <- permute::shuffleSet(ncol, nset=nrow, quiet=TRUE);
			perm <- t(sapply(seq_len(nrow(P)), function(i, P, M) M[i,P[i,]],M=M,P=P))
			D <- c(D,as.vector(distfun(t(perm))))
		}
		threshold <- quantile(D,probs=0.05)
	
		my_dists <- distfun(t(my_profiles))
#		my_dist_signif <- as.matrix(my_dists) < (quantile(D, probs=0.05/prod(dim(my_dists))/2))
		pairwise_threshold <- quantile(D, probs=0.05/prod(dim(my_dists))/2)
		my_hclust <- hclustfun(my_dists)
		my_sig <- cutree(my_hclust, h=threshold)
		perm_signif <- rainbow(n=max(my_sig))[my_sig]
	}
	# Add pair-wise cluster significances.

	## Set-up Heatmap ##
	# heatmap colours & binning
	heatcols <- colorRampPalette(c("blue","white","red"))(255)
	#heatcols <- colorRampPalette(c("white","black"))(255)

	# Column bar 1 = matches
	matched <- factor_counts(matches$labels);
	unmatched <- names(matched[matched == 1]);
	match_lab <- as.character(matches$labels)
	match_lab[match_lab %in% unmatched] <- "00_Unmatched"
	match_lab <- factor(match_lab)
	match_col = c("white",RColorBrewer::brewer.pal(length(levels(match_lab))-1, "Set2")) # Change "white" to background color
	Matches <- match_col[match_lab];
	ColumnCols = data.frame(Matched=Matches)
	# Column bar 2 = permute threshold;
	if (npermute > 0) {
		ColumnCols$Signif <- perm_signif
	}
	# Column bar 3 = ann
	if (!is.null(ann)) {
		ann <- as.factor(ann)
		#ann_col <- RColorBrewer::brewer.pal(length(levels(ann)), "Set3")
		ann_col <- colorRampPalette(c("#d9d9d9","#fccde5","#bebada","#bc80bd","#80b1d3","#8dd3c7","#b3de69","#ccebc5","#ffffb3","#ffed6f","#fdb462","#fb8072"))(length(levels(ann)))
		ColumnCols$Known <- ann_col[ann];
	}

	if (show_genes) {
		heatout <- heatmap.3(my_profiles, trace="n", scale="row", col=heatcols, symbreaks=TRUE, key.title="", key.xlab="Relative Expression", hclustfun=hclustfun, distfun=distfun, ColSideColors=as.matrix(ColumnCols), ColSideColorsSize=length(ColumnCols[1,]))
	} else {
		pro_vs_pro = distfun(t(my_profiles))
		if (exists("threshold")) {
			my_max = max(pro_vs_pro);
			tmp = as.matrix(pro_vs_pro)
			tmp[tmp > threshold] = my_max
			pro_vs_pro = tmp;
		}
		dendro <- hclustfun(as.dist(pro_vs_pro))

		heatout <- heatmap.3(pro_vs_pro, trace="n", scale="none", col=rev(heatcols), symbreaks=FALSE,  key.title="", key.xlab="Distance", ColSideColors=as.matrix(ColumnCols), ColSideColorsSize=length(ColumnCols[1,]), symm=TRUE, Rowv=as.dendrogram(dendro), Colv=as.dendrogram(dendro))
	}

	return(invisible(list(heatmap_out=heatout, sig_groups = perm_signif, dist_mat=pro_vs_pro, sig_threshold=threshold)))
}
