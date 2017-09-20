# Error Catching needs: duplicate row names, no consensus genes.
# Error Catching needs: NAs in profiles
# Error Catching needs: cell_type levels with no samples.


# Generalizable - Step 1 = Get profiles
get_cluster_profile <- function(expr_mat, clusters, norm_method="none", is.log=FALSE, out.log=2, feature_selection=marker.features) {
	# Add quantile normalization somewhere????
	if (min(factor_counts(clusters)) < 2) {stop("Error: Cannot have a cell-type with a single sample.")}
	# is.log should hold the base of the log if it has been log-transformed.
	# If no norm here, need at least FS!
	if (is.log) {
		expr_mat <- is.log^expr_mat;
		expr_mat <- expr_mat - 1;
	}
	

	profiles <- my_row_mean_aggregate(expr_mat, clusters)
	is.feature <- feature_selection(expr_mat, clusters) # now matrix of same dimensions of profiles.
	if (norm_method=="CPM") {
		# CPM
		normfact <- colSums(profiles);
	} else if (norm_method %in% c("TMM", "RLE", "upperquartile", "none")) {
		# TMM
		#require("edgeR")
		normfact <- edgeR::calcNormFactors(as.matrix(profiles), method=norm_method)
	} else if (norm_method == "quantile") {
		profiles <- preprocessCore::normalize.quantiles(profiles)
		expr_mat <- preprocessCore::normalize.quantiles(expr_mat, copy=FALSE)
	} else {
		warning("Warning: Unrecognized normalization method. Profiles will be scaled but not normalized.")
		normfact <- rep(1, times=length(profiles[1,]))
	}

	# Apply & scale to ~ 1,000,000 counts per profile
	profiles <- t(t(profiles)/normfact)
	scalefact <- 1000000/median(colSums(profiles))
	profiles <- profiles*scalefact

	expr_mat_scalefact <- 1000000/median(colSums(expr_mat))
	expr_mat <- expr_mat*expr_mat_scalefact

	 # out.log should hold the base of the log if should be log-transformed.
	if (out.log) {
		profiles <- log(profiles+1)/log(out.log);
		expr_mat <- log(expr_mat+1)/log(out.log);
	}

	return(list(norm_mat = expr_mat, profiles=profiles, is.feature = is.feature))
}

# Generalizable - Step 2 = Match clusters using reciprocal best-hits
get_reciprocal_hits <- function(profiles1, profiles2, features1, features2) {
	# features = gene names of features to use

	g1 <- rownames(profiles1)
	g2 <- rownames(profiles2)

        if ( (length(g1) != length(g2)) | (suppressWarnings(sum(g1 == g2)) != length(g1)) ) {
		keep_g <- rownames(profiles1) %in% rownames(profiles2)
                profiles1 <- profiles1[keep_g,]
		features1 <- features1[keep_g,]
		g_match <- match(rownames(profiles1), rownames(profiles2))
                profiles2 <- profiles2[g_match,]
		features2 <- features2[g_match,]
#                dim1 <- dim(profiles1);
#                dim2 <- dim(profiles2);
        }
#        features1 <- profiles1[,dim1[2]]
#        features2 <- profiles2[,dim2[2]]
#        profiles1 <- as.matrix(profiles1[,-dim1[2]])
#        profiles2 <- as.matrix(profiles2[,-dim2[2]])

#      	if (is.na(features)) { # throw warning when features is not NA
#      	        # feature of one or other dataset.
#      	        consensus_features = features1==1 | features2==1;
#      	} else {
#      	        # user provided features.
#      	        consensus_features = rownames(profiles1) %in% features;
#      	}

	matches1_IDs <- rep(-1, times=ncol(profiles1))
	matches1_vals <- rep(-1, times=ncol(profiles1))
	matches2_IDs <- rep(-1, times=ncol(profiles2))
	matches2_vals <- rep(-1, times=ncol(profiles2))

	for (i in 1:ncol(profiles1)) {
		for(j in 1:ncol(profiles2)) {
			# feature olap
			cons_features <- features1[,i] | features2[,j]
			sim_score <- cor(profiles1[cons_features,i], profiles2[cons_features, j], 
					method="spearman")
			if (sim_score > matches1_vals[i]) {
				matches1_vals[i] = sim_score
				matches1_IDs[i] = j
			}
			if (sim_score > matches2_vals[j]) {
				matches2_vals[j] = sim_score
				matches2_IDs[j] = j
			}
		}
	}

        # Get Matches
#        match_cors <- cor(profiles1[consensus_features,], profiles2[consensus_features,], method="spearman");
#        matches2 = apply(match_cors, 2, function(x){which(x==max(x))})

        matches1 = matches1_IDs
        matches2 = matches2_IDs

        # Get Reciprocal Best-Hits
        recip = which(1:length(matches1) == matches2[matches1])
        recip_table = cbind(recip,matches1[recip])
        colnames(recip_table) = c("profiles1", "profiles2")
#        match_table = cbind(matches1, matches2) # not necessarily of equal length
#        colnames(match_table) = c("profiles1->profiles2", "profiles2->profiles1")
        return(list(matches1=matches1, matches2=matches2, recip=recip_table))
}

get_multi_hits <- function (profiles1, profiles2, features1, features2, sig.threshold=0.05, CI.level=0.95) {
        # features = gene names of features to use
	g1 <- rownames(profiles1)
	g2 <- rownames(profiles2)

        if ( length(g1) != length(g2) | suppressWarnings(sum(g1 == g2)) != length(g1) ) {
                keep_g <- rownames(profiles1) %in% rownames(profiles2)
                profiles1 <- profiles1[keep_g,]
                features1 <- features1[keep_g,]
                g_match <- match(rownames(profiles1), rownames(profiles2))
                profiles2 <- profiles2[g_match,]
                features2 <- features2[g_match,]
        }

        cor_matrix <- matrix(-1, ncol=ncol(profiles2), nrow=ncol(profiles1))
        cor_CI_high_matrix <- matrix(-1, ncol=ncol(profiles2), nrow=ncol(profiles1))
        for (i in 1:ncol(profiles1)) {
                for(j in 1:ncol(profiles2)) {
                        # feature olap
                        olap <- sum(features1[,i] & features2[,j])
                        pval <- phyper(olap, sum(features1[,i]), sum(features1[,i]==0),
                                        sum(features2[,j]), lower.tail=FALSE)
                        if(pval > sig.threshold){next;}
                        cons_features <- features1[,i] | features2[,j]
                        sim_score <- suppressWarnings(cor.test(profiles1[cons_features,i], profiles2[cons_features, j],
                                        method="spearman")) # need to suppress warnings.
                        if (sim_score$p.value > sig.threshold) {next;}
                        cor_matrix[i,j] <- sim_score$estimate
                        # CI based on Fieller Hartley Pearson 1957 & confirmed by Bishara & Hittner 2017
                        # equations from Ruscio 2008
                        zr <- atanh(sim_score$estimate)
                        CI_hi <- tanh(zr+sqrt(1/(sum(cons_features)-3)*qnorm((1+CI.level)/2)))
                        cor_CI_high_matrix[i,j] <- CI_hi
                }
        }
	rownames(cor_matrix) <- colnames(profiles1)
	colnames(cor_matrix) <- colnames(profiles2)
	top_hits1 <- apply(cor_matrix, 1, max)
	all_hits1 <- cor_CI_high_matrix > top_hits1
	hit_tab1 <- which(all_hits1, arr.ind=TRUE)
	top_hits2 <- apply(cor_matrix, 2, max)
	all_hits2 <- cor_CI_high_matrix > top_hits2
	hit_tab2 <- which(all_hits2, arr.ind=TRUE)

	recip_table <- rbind(hit_tab1, hit_tab2)
        colnames(recip_table) = c("profiles1", "profiles2")
	recip_table <- recip_table[duplicated(recip_table),] # require reciprocal

	return(list(matches1=NA, matches2=NA, recip=recip_table))
}


combine_and_match_clusters <- function(profile_list, features=NA, multihit=FALSE, sig.threshold=0.05, CI.level=0.95, suppress.plot=TRUE){
	# features = gene names of features to use
	if (length(profile_list) < 2) {stop("Error: Cannot combine & match profiles from a single dataset")}

	# Dataset names
	if (is.null(names(profile_list))) {
		warning("Warning: no dataset names, creating default names");
		names(profile_list) <- paste("Data", 1:length(profile_list));
	}
	# attach dataset names to cluster names
	for (i in 1:length(profile_list)) {
		colnames(profile_list[[i]]$profiles) <- paste(names(profile_list)[i], colnames(profile_list[[i]]$profiles), sep="_");
		reorder = order(rownames(profile_list[[i]]$profiles));
		profile_list[[i]]$profiles <-  profile_list[[i]]$profiles[reorder,]
		profile_list[[i]]$norm_mat <-  profile_list[[i]]$norm_mat[reorder,]
	}
	
	# for each pair of datasets get reciprocal best-hits
	pairs = combn(names(profile_list), 2)
	RECIP <- vector()
	for(j in 1:length(pairs[1,])) {
		profile_obj1 = profile_list[[pairs[1,j]]]
		profile_obj2 = profile_list[[pairs[2,j]]]
		features1 = profile_obj1$is.feature
		features2 = profile_obj2$is.feature
		if (!is.na(features)) {
			features1 = matrix(rep(features, times=ncol(features1)), ncol=ncol(features1));
			features2 = matrix(rep(features, times=ncol(features2)), ncol=ncol(features2));
		}

		if (multihit) {
			out = get_multi_hits(profile_obj1$profiles, profile_obj2$profiles, 
					features1=features1, features2=features2,
					sig.threshold=sig.threshold, CI.level=CI.level)
		} else {
			out = get_reciprocal_hits(profile_obj1$profiles, profile_obj2$profiles, 
					features1=features1, features2=features2)
		}
		if (is.null(dim(out$recip))) {
			recip = as.data.frame(t(out$recip))
		} else {
			recip = as.data.frame(out$recip)
		}
		recip[,1] <- colnames(profile_list[[pairs[1,j]]]$profiles)[recip[,1]]
		recip[,2] <- colnames(profile_list[[pairs[2,j]]]$profiles)[recip[,2]]
		RECIP <- rbind(RECIP, recip)
	}
	# match them all using graph connected components
	#require("igraph")
	my_graph <- igraph::graph_from_edgelist(as.matrix(RECIP), directed=FALSE)
	matched_clusters <- igraph::components(my_graph)$membership
	if (!suppress.plot) {
		plot(my_graph) # colour this by dataset
	}
	
	# Combine profiles together
	# consensus genes
	consensus_genes = rownames(profile_list[[1]]$profiles)
	for (i in 2:length(profile_list)) {
		consensus_genes <- consensus_genes[consensus_genes %in% rownames(profile_list[[i]]$profiles)]	
	}
	# combine profiles
# Need to check & throw out datasets which have no matches with the other datasets.
	# duplicate code
	profile_tab1 <- profile_list[[1]]$profiles
	these_profiles <- as.matrix(profile_tab1[ rownames(profile_tab1) %in% consensus_genes, ])
	colnames(these_profiles) <- colnames(profile_tab1)
	Combined <- these_profiles
	these_features <- profile_list[[1]]$is.feature[rownames(profile_tab1) %in% consensus_genes, ]
	colnames(these_features) <- colnames(profile_tab1)
	Combined_features <- these_features
	# end duplicate code
	for (i in 2:length(profile_list)) {
		# duplicate code
		profile_tabi <- profile_list[[i]]$profiles
		these_profiles <- as.matrix(profile_tabi[rownames(profile_tabi) %in% consensus_genes, ])
		colnames(these_profiles) <- colnames(profile_tabi)
		Combined <- cbind(Combined, these_profiles)
		these_features <- profile_list[[i]]$is.feature[rownames(profile_tabi) %in% consensus_genes, ]
		colnames(these_features) <- colnames(profile_tabi)
		Combined_features <- cbind(Combined_features, these_features)
		#end duplicate code
	}
	#colnames(Combined_features) <- names(profile_list)

	# merge match clusters
	Combined_Clusters <- colnames(Combined)
	for (comp in 1:max(matched_clusters)) {
		Combined_Clusters[Combined_Clusters %in% names(matched_clusters)[matched_clusters==comp]] <- paste("Match",comp, sep="_")
	}
	# dataset labels
	Combined_Dataset <- vector();
	for (i in 1:length(profile_list)) {
		# Check if all celltypes for this dataset are present in "Combined Clusters" if so then no cross-dataset matches & dataset should be excluded.
		profile_tab <- profile_list[[i]]$profiles
		these_groups <- colnames(profile_tab)
		if (sum(these_groups %in% Combined_Clusters) == length(these_groups) ) {
			# No cross-dataset matches
			warning(paste("Datasets: ", names(profile_list)[i], "shares no cell-type with other datasets and will be excluded.", sep=""))
			exclude <- Combined_Clusters %in% these_groups;
			Combined_Clusters <- Combined_Clusters[!exclude]
			Combined <- Combined[,!exclude]
			Combined_features <- Combined_features[!exclude]
		} else {
			labs <- rep(names(profile_list)[i], times=ncol(profile_tab))
			Combined_Dataset <- c(Combined_Dataset, labs)
		}
	}
	if (length(unique(Combined_Dataset)) == 1) {
		warning("Warning: Only one dataset remains after combining!\nYou may want to try again with a higher value for sig.threshold.")
	}

	return(list(profiles = Combined, labels = Combined_Clusters, dataset = Combined_Dataset, features = Combined_features, match_graph = my_graph))
}

# Generalizable - Step 3 = GLM to calculate batch effects 
#		& Step 4 = remove batch effects from all profiles

glm_of_matches <- function(matches) { #This is not optimized
	# This should work for splitting cell-type specific DE & batch effect DE
	# But I need return more than just the effect sizes
	matched <- factor_counts(matches$labels);
	matched <- names(matched[matched >1])

	matched_matrix <- matches$profiles[,matches$labels %in% matched]
	matched_type <- matches$labels[matches$labels %in% matched]
	matched_batch <- matches$dataset[matches$labels %in% matched]

	batch_lvl = levels(factor(matches$dataset))
	type_lvl = levels(factor(matched_type))
	all_matrix <- matches$profiles
	# I think this former case is mathematically impossible.
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

		# rowID based on batchID for glm output
		# 1 = intercept, 2-x = cell-types(-1), x-n = batches
		effect_vec = rep(1, times=length(matched_batch)); # matched only
		all_effect <- rep(1, times=length(matches$dataset)) 
	
		for (i in 2:length(batch_lvl)) {
			effect_vec[matched_batch == batch_lvl[i]] <- i+length(type_lvl);
			all_effect[matches$dataset == batch_lvl[i]] <- i+length(type_lvl);
		}	
		calc_batch_effect <- function(x) {
			out <- glm(x~matched_type+matched_batch); # Do I force the intercept to be zero i.e. ~0+matched_type+matched_batch
			sout <- summary(out)$coeff;
			#effects <- c(0,sout[,1]);
			#pval <- c(0,sout[,4]);
			##effects = c(0,glm(x~matched_type+matched_batch)$coef); # Do I force the intercept to be zero i.e. ~0+matched_type+matched_batch
			#return(rbind(effects, pval)) # this no longer make a table
			return(sout)
		}
		batch_effects <- apply(matched_matrix, 1, calc_batch_effect)
		be_dims <- dim(batch_effects)
		effect_size <- batch_effects[1:(be_dims[1]/4),] # top quarter
		effect_size <- rbind(rep(0, times=ncol(effect_size)), effect_size)
		p_values <- batch_effects[be_dims[1]-(rev(1:(be_dims[1]/4))-1),] # bottom quarter
		p_values <- rbind(rep(0, times=ncol(p_values)), p_values)

		effect_names <- c("Reference", "Intercept", levels(factor(matched_type))[-1], levels(factor(matched_batch))[-1])
		rownames(effect_size) <- effect_names
		rownames(p_values) <- effect_names

		corrected <- all_matrix - t(effect_size[all_effect,])
	}
	return(list(corrected_profiles=corrected, batch_effects=effect_size, effect_p_values=p_values))
}

# Generalizable - Weighted Means

wmeans_of_matches <- function (matches){ #This is not optimized
	all_matrix <- matches$profiles
	dataset_fac <- factor(matches$dataset)
	batch_effects <- my_row_mean_aggregate(all_matrix, dataset_fac);
	
	corrected <- all_matrix-batch_effects[,dataset_fac]
	
	return(list(corrected_profiles=corrected, batch_effects=batch_effects))
}

# Generalizable - Step 5 = Cluster corrected profiles
#		& Step 6 = Calculate meta-profiles
cluster_profile_heatmap <- function(corrected_profiles, matches, features_only=TRUE, show_genes=FALSE, npermute=0, distfun=function(x){as.dist(1-cor(t(x), method="spearman"))}, hclustfun=function(x){hclust(x,method="complete")}, ann=NULL, dataset=NULL) {
	#source("/nfs/users/nfs_t/ta6/R-Scripts/heatmap.3.R")

	my_profiles <- corrected_profiles
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
	perm_signif <- NA;
	if (npermute > 0) {
		set.seed(101); # reproduciblity
		for(rep in 1:npermute) {
			P <- permute::shuffleSet(ncol, nset=nrow, quiet=TRUE);
			perm <- t(sapply(base::seq_len(nrow(P)), function(i, P, M) M[i,P[i,]],M=M,P=P))
			D <- c(D,as.vector(distfun(t(perm))))
		}
		threshold <- quantile(D,probs=0.05)
	
		my_dists <- distfun(t(my_profiles))
#		my_dist_signif <- as.matrix(my_dists) < (quantile(D, probs=0.05/prod(dim(my_dists))/2))
		pairwise_threshold <- quantile(D, probs=0.05/prod(dim(my_dists))/2)
		my_hclust <- hclustfun(my_dists)
		my_sig <- cutree(my_hclust, h=threshold)
		perm_signif <- rainbow(n=max(my_sig))[my_sig]
	} else {
		threshold=NA
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
	match_col = c("white", colorRampPalette(c("#ffffcc","#fed9a6","#fbb4ae","#decbe4","#b3cde3","#ccebc5"))(length(levels(match_lab))-1))
#	match_col = c("white",RColorBrewer::brewer.pal(length(levels(match_lab))-1, "Set2")) # Change "white" to background color
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
	if (!is.null(dataset)) {
		dataset <- as.factor(dataset)
		#ann_col <- RColorBrewer::brewer.pal(length(levels(ann)), "Set3")
		data_col <- colorRampPalette(c("#fdb462", "#e5c494","#b3de69","#ffd92f","#fccde5"))(length(levels(dataset)))
		ColumnCols$Dataset <- data_col[dataset];
	}

	if (show_genes) {
		heatout <- heatmap.3(my_profiles, trace="n", scale="row", col=heatcols, symbreaks=TRUE, key.title="", key.xlab="Relative Expression", hclustfun=hclustfun, distfun=distfun, ColSideColors=as.matrix(ColumnCols), ColSideColorsSize=length(ColumnCols[1,]))
	} else {
		pro_vs_pro = distfun(t(my_profiles))
		if (exists("threshold") & !is.na(threshold)) {
			my_max = max(pro_vs_pro);
			tmp = as.matrix(pro_vs_pro)
			tmp[tmp > pairwise_threshold] = my_max
			tmp <- tmp/my_max; # rescale to 0-1.
			pro_vs_pro = tmp;
		} else {
			threshold <- NA
			my_max = max(pro_vs_pro);
			tmp = as.matrix(pro_vs_pro)
			tmp <- tmp/my_max; # rescale to 0-1.
			pro_vs_pro = tmp;
		}
		dendro <- hclustfun(as.dist(pro_vs_pro))

		heatout <- heatmap.3(pro_vs_pro, trace="n", scale="none", col=rev(heatcols), symbreaks=FALSE,  key.title="", key.xlab="Distance", ColSideColors=as.matrix(ColumnCols), ColSideColorsSize=length(ColumnCols[1,]), symm=TRUE, Rowv=as.dendrogram(dendro), Colv=as.dendrogram(dendro))
	}

	return(invisible(list(heatmap_out=heatout, sig_groups = perm_signif, dist_mat=pro_vs_pro, sig_threshold=threshold)))
}

# Use the norm matrix from the new profiles method.
correct_sng_cells <- function(norm_mat, dataset_name, glm_out, allow.negatives=FALSE) { # SLOW!
	#dataset_row = grep(paste("matched_batch",dataset_name, sep=""), rownames(glm_out$batch_effects))
	dataset_row = grep(dataset_name, rownames(glm_out$batch_effects))
	if (length(dataset_row)==0) {
		print("This is the reference batch - no correction to be made")
		return(norm_mat);
	} 
	keep_genes <- colnames(glm_out$batch_effects)
	norm_mat <- norm_mat[rownames(norm_mat) %in% keep_genes,]

	if (allow.negatives) {
		correct <- t(glm_out$batch_effects[rep(dataset_row, times=ncol(norm_mat)),])
		return(norm_mat-correct);
	} else {
		change = glm_out$batch_effects[dataset_row,]
		increases = change < 0
		norm_mat[increases,] <- norm_mat[increases,]-change[increases]
		
		for(i in which(change > 0)) {
			mean_change = change[i]
			total_change = mean_change*ncol(norm_mat);
			if (total_change >= sum(norm_mat[i,])) {
				# Short circuit if total change is greater than 
				# turning whole row to zero.
				norm_mat[i,] <- rep(0, length=ncol(norm_mat));
				next;
			}

			# This is trick b/c [total_change-sum(expr_mat[i,expr_mat[i,] < x)]/sum(expr_mat[i,] >=x) < min(expr_mat[i,] >=x);
			threshold = norm_mat[i,]
			zeros_change <- sapply(threshold, function(x) {sum(norm_mat[i,norm_mat[i,] < x])})
			new_mean_change <- (total_change-zeros_change)/sapply(threshold, function(y){sum(norm_mat[i,] >= y)})
			threshold_correct = min(threshold[new_mean_change < threshold])

			zeros_change <- sum(norm_mat[i,norm_mat[i,] < threshold_correct])
			norm_mat[i,norm_mat[i,] < threshold_correct] <- 0;

			remaining_change <- total_change-zeros_change
		
			norm_mat[i,] <- norm_mat[i,]-remaining_change/sum(norm_mat[i,] > 0)
			norm_mat[i,norm_mat[i,] < 0 ] <- 0
		}
		
	}
	return(norm_mat);
}
