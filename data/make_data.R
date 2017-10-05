require("RColorBrewer")
require("gplots")
require("scater")
require("CellTypeProfiles")

map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Mmus_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}
get_ICM_TE_assignments<- function(x, labels) {
        data <- log(x+1)/log(2);
        # TE & ICM assignment
        TE_markers=c("Elf5","Eomes","Cdx2")
        ICM_markers = c("Sox2","Pou5f1","Nanog")
        scale <- function(x) {(x-mean(x))/sd(x)}
        ICMscore = rowMeans(apply(data[rownames(data) %in% ICM_markers,], 1, scale))
        TEscore = rowMeans(apply(data[rownames(data) %in% TE_markers,], 1, scale))
        blasts = (grepl("blast", labels) | grepl("32cel", labels))
        new_lab = as.character(labels);
        if (sum(blasts) > 0) {
                new_lab[blasts & ICMscore > TEscore] = "ICM"
                new_lab[blasts & ICMscore < TEscore] = "TE"
        }
        names(new_lab) <- colnames(x)
        return(new_lab);
}

human_ICM_TE_assignments <- function(data, labels) {
        # TE & ICM assignment
        TE_markers=c("CDX2", "POU5F1", "GATA2", "GATA3", "CLDN10")
	EPI_markers = c("SOX2", "TDGF1", "DPPA5", "GDF3", "PRDM14")
	PE_markers = c("PDGFRA", "HNF1B", "BMP2", "FGFR2", "KIT")
        scale <- function(x) {(x-mean(x))/sd(x)}
        PEscore = rowMeans(apply(data[rownames(data) %in% PE_markers,], 1, scale))
        EPIscore = rowMeans(apply(data[rownames(data) %in% EPI_markers,], 1, scale))
        TEscore = rowMeans(apply(data[rownames(data) %in% TE_markers,], 1, scale))
        blasts = (grepl("blast", labels) | grepl("32cel", labels))
        new_lab = as.character(labels);
        if (sum(blasts) > 0) {
                new_lab[blasts & PEscore > EPIscore & PEscore > TEscore] = "PE"
                new_lab[blasts & EPIscore > PEscore & EPIscore > TEscore] = "EPI"
                new_lab[blasts & TEscore > PEscore & TEscore > EPIscore] = "TE"
        }
        names(new_lab) <- colnames(data)
        return(new_lab);
}
convert_to_integer <- function(mat) {
        mat <- round(as.matrix(mat))
        storage.mode(mat) <- "integer"
        mat = mat[rowSums(mat) > 0,]
        return(mat)
}
get_log_norm <- function(SCE) {
        dat <- counts(SCE);
        fac <- colSums(dat);
        norm <- t(t(dat)/fac*median(fac))
        return(log(norm +1)/log(2));
}



source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/Consistent_Setup.R")

# Deng
load("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Deng_embryo_clean.RData")
deng_count = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = FALSE)
deng_list = normalize_data(Deng_embyro_list$data, labels = Deng_embyro_list$labels, is.counts = TRUE)
deng_lab <- get_ICM_TE_assignments(deng_list$data, Deng_embyro_list$labels)
deng_lab[deng_lab=="early2cell"] <- "zygote"
deng_lab[deng_lab=="mid2cell"] <- "2cell"
deng_lab[deng_lab=="late2cell"] <- "2cell"
Deng_profile = get_cluster_profile(deng_list$data, deng_lab, norm_method="CPM", is.log=FALSE, out.log=2, feature_selection=marker.features)
Deng_profile$norm_mat <- NULL

# Zhong
zhong = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE57249_fpkm_ZHONG.txt", header=TRUE);
zhong = zhong[!duplicated(zhong[,1]),]
rownames(zhong) = zhong[,1]
zhong = zhong[,2:length(zhong[1,])]
zhong = as.matrix(zhong);

zhong_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE57249_labels_ZHONG.txt");
ultralow = which(rowMeans(zhong) < 10^-5)
zhong = zhong[-ultralow,]
zhong_list = normalize_data(zhong, labels=zhong_labels, is.counts=FALSE)
zhong_count = convert_to_integer(zhong_list$data);
zhong_list$data <- zhong_list$data[rownames(zhong_list$data) %in% rownames(zhong_count),]
#zhong_lab <- get_ICM_TE_assignments(zhong_list$data, unlist(zhong_labels))
zhong_lab=as.character(unlist(zhong_labels))

Zhong_profile = get_cluster_profile(zhong_list$data, zhong_lab, norm_method="CPM", is.log=FALSE, out.log=2, feature_selection=marker.features)
zhong_list$labels <- zhong_lab;
biase_list <- zhong_list
save(biase_list, file="Biase.rda")

# Xue
Xue_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE44183_mouse_expression_mat.txt", header=TRUE)
Xue_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE44183_mouse_expression_labels.txt", header=FALSE)
Xue_labels = as.character(unlist(Xue_labels))
Xue_labels[Xue_labels == "2cellmixed"] = "2cell"
Xue_labels[Xue_labels == "4cellmixed"] = "4cell"
Xue_labels[Xue_labels == "8cellmixed"] = "8cell"
Xue_labels[Xue_labels == "Pronucleus"] = "zygote"
Xue_labels[Xue_labels == "Oocyte"] = "oocyte"
Xue_labels[Xue_labels == "Morula"] = "16cell"
Xue_list = normalize_data(Xue_data, labels=Xue_labels, is.counts=FALSE)
Xue_count = convert_to_integer(Xue_list$data);
Xue_list$data <- Xue_list$data[rownames(Xue_list$data) %in% rownames(Xue_count),]
Xue_lab <- get_ICM_TE_assignments(Xue_list$data, Xue_labels)

Xue_profile = get_cluster_profile(as.matrix(Xue_list$data), Xue_lab, norm_method="CPM", is.log=FALSE, out.log=2, feature_selection=marker.features)
Xue_profile$norm_mat <- NULL

# Fan
Fan_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE53386_matrix_fpkms.tsv", header=TRUE)
Fan_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/GSE53386_labels.tsv", header=FALSE)
Fan_labels = as.character(unlist(Fan_labels))
Fan_labels[Fan_labels == "a-AM_treated_2-cell"] = "oocyte"
Fan_labels[Fan_labels == "2-cell"] = "2cell"
Fan_labels[Fan_labels == "4-cell"] = "4cell"
Fan_labels[Fan_labels == "8-cell"] = "8cell"
Fan_labels[Fan_labels == "morula"] = "16cell"
Fan_list = normalize_data(Fan_data, labels=Fan_labels, is.counts=FALSE)
Fan_count = convert_to_integer(Fan_list$data)
Fan_list$data <- Fan_list$data[rownames(Fan_list$data) %in% rownames(Fan_count),]
Fan_lab <- get_ICM_TE_assignments(Fan_list$data, Fan_labels)
Fan_profile = get_cluster_profile(as.matrix(Fan_list$data), Fan_lab, norm_method="CPM", is.log=FALSE, out.log=2, feature_selection=marker.features)
Fan_profile$norm_mat <- NULL

# Goolam
Goolam_data = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Goolam_et_al_2015_count_table.tsv", header=T)
Goolam_labels = read.table("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Goolam_et_al_2015_count_table_labels.tsv", header=F)
Goolam_list = normalize_data(Goolam_data, unlist(Goolam_labels[1,]), is.counts=TRUE)
Goolam_counts = normalize_data(Goolam_data, unlist(Goolam_labels[1,]), is.counts=FALSE)
Goolam_list$labels = (as.character(Goolam_list$labels))
Goolam_list$labels[Goolam_list$labels =="32cell"] = "blast"
rownames(Goolam_counts$data) = ensg2symbol(rownames(Goolam_counts$data))
rownames(Goolam_list$data) = ensg2symbol(rownames(Goolam_list$data))
Goolam_lab <- get_ICM_TE_assignments(Goolam_list$data, Goolam_list$labels)
Goolam_profile = get_cluster_profile(Goolam_list$data, Goolam_lab, norm_method="CPM", is.log=FALSE, out.log=2, feature_selection=marker.features)
Goolam_profile$norm_mat <- NULL

Biase_profile <- Zhong_profile;

truth <- c(levels(factor(deng_lab)), levels(factor(zhong_lab)), levels(factor(Xue_lab)), levels(factor(Fan_lab)), levels(factor(Goolam_lab)))
save(Deng_profile, Biase_profile, Xue_profile, Fan_profile, Goolam_profile, truth, file="devo_profiles.rda")
