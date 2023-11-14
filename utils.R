dpe_enrich_plot  <- function(m, M1, myVar, pthres = 0.01, addLabel = F, leadThres = 5) {
	#double check factor levels
	contr.string <- paste0(myVar, " - DiagnosisControl")
	#print(contr.string)
levels(as.factor(meta_$Diagnosis))
#Control is reference
fitM1 <- limma::lmFit(log1p(t(m)), M1)
#contrast <- do.call(limma::makeContrasts, list(PatientVScontrol = contr.string, levels = M1))
#colnames(contrast) <- "PatientVScontrol"
#print(contrast)
#fitM1 <- limma::contrasts.fit(fitM1, contrast) 
fitM1  <- limma::eBayes(fitM1, trend = T, robust = T)
resM1 <- limma::topTable(fitM1, myVar, n = Inf) %>%
	#limma::topTable(fitM1, "PatientVScontrol", n = Inf) %>%
	as.data.frame %>%
	mutate(Accession = rownames(.)) %>%
	left_join(., 
		  dataDF %>% select(Accession, gene_name, gene_id),
		  by = "Accession")
print("Up and down:")
print(resM1 %>% filter(adj.P.Val < 0.05) %>% 
	group_by(sign(logFC)) %>%
	tally)
resM1 %>% filter(adj.P.Val < 0.05) %>% 
 	ggplot(., aes( x = logFC)) +
	geom_histogram() +
	theme_cowplot() -> p1
resM1 %>%  
 	ggplot(., aes( x = P.Value)) +
	geom_histogram() +
	theme_cowplot() -> p2

resM1 <- resM1 %>%
#	filter(adj.P.Val < 0.05) %>%
	arrange(desc(abs(logFC))) %>% 
	mutate(Rank = seq(1, nrow(.))) %>%
	mutate(gene_name = ifelse(Accession == "P02462", "COL4A1", gene_name)) %>%
	mutate(Type = case_when(
				gene_name %in% oxphos$gene_name ~ "OXPHOS",
				grepl("^NDUFS", gene_name) ~ "OXPHOS",
				grepl("^SLC", gene_name) ~ "Solute carrier",
				grepl("^PSM", gene_name) ~ "Proteasome",
				grepl("^SNC", gene_name) ~ "Synuclein",
				grepl("COL", gene_name) ~ "Collagens",
				grepl("^MR", gene_name) | gene_name == "DAP3" ~ "Mitochondrial ribosome",
				#gene_name == "DAP3" ~ "Mitochondrial ribosome",
				grepl("RPL", gene_name) ~ "Ribosomal",
				grepl("RPS", gene_name) ~ "Ribosomal",
				grepl("HSP", gene_name) ~ "Heatshock protein",
				grepl("ATP", gene_name) ~ "ATP",
				grepl("SNCA", gene_name) ~ "Alpha-syn",
				grepl("KRT", gene_name) ~ "Keratine gene family",
				grepl("NAD", gene_name) ~ "NAD in gene name",
				gene_name %in% c("MARK1", "MAPK3", "MAP2K6") ~ "Kinase",
				gene_name %in% (mgp %>%
						filter(cell_type == "Neurons") %>%
						pull(gene_name)) ~ "Neuron marker",
				gene_name %in% (mgp %>%
						filter(cell_type == "Oligodendrocytes") %>%
						pull(gene_name)) ~ "Oligo marker",
				gene_name %in% (mgp %>%
						filter(cell_type == "Astrocytes") %>%
						pull(gene_name)) ~ "Astrocyte marker",
				gene_name %in% (mgp %>%
						filter(cell_type == "Microglia") %>%
						pull(gene_name)) ~ "Microglia marker",
				gene_name %in% synapses ~ "Synapse",
				grepl("SYT", gene_name) ~ "Synapse",
				grepl("^CCT", gene_name) ~ "Chaperonin containing CCT",
				TRUE ~ "Other")) 
resM1 %>%
	mutate(Type = as.factor(Type)) %>%
	dplyr::group_by(Type) %>%
	dplyr::summarize(medianRank = median(Rank)) %>%
	arrange(medianRank) %>%
	pull(Type) -> orderType
resM1 %>%
	mutate(Type = factor(Type, levels = unique(orderType))) %>%
	mutate(label = ifelse((Type == "Other" &  Rank < 10) | gene_name %in% haris_g, gene_name, "")) %>%
	filter(adj.P.Val < 0.05) %>%
	ggplot(., aes(x = Type, y = Rank, col = Type)) +
	ggbeeswarm::geom_beeswarm(size =.6, dodge.width = 0.9 ) +
	coord_flip() +
	geom_hline(yintercept = 50, lty = "dashed", col = "grey") +
	geom_boxplot(fill = "transparent") +
	theme_cowplot() -> p3

if (addLabel) {
	p3 <- p3 + 
		ggrepel::geom_label_repel(aes(label = label), max.overlaps = Inf)
}
gp <- (p1 + p2) / p3
#Enrichment
resM1 %>%   
dplyr::arrange(desc(t)) %>%
	mutate(score = -log10(P.Value)*sign(logFC)) %>%
	dplyr::arrange(desc(abs(score))) %>%

    # Im doing these in case there is duplicates in gene names
    # This might remove genes you are actually interested in due to them having two unique uniprot ids (Accession) mapped to the same gene_name (could be isoforms) So Id recommend looking at duplicates and see if some of them are of interest to you
    dplyr::distinct(gene_name, .keep_all=TRUE) %>%
    # Same goes for this
    dplyr::filter(!is.na(gene_name)) %>%
    as.data.frame -> res_df
rownames(res_df) <- res_df$gene_name
ranks <- setNames(res_df$t, res_df$gene_name)
ranks <- setNames(res_df$score, res_df$gene_name)

gseas.PD <- lapply(names(pathways), function(pathwayName){
    message(paste0("Running ", pathwayName, "..."))
    fgsea(pathways[[pathwayName]],
	  ranks,
	  minSize = 20,
	  maxSize = 200,
	  eps = 0) %>%
        as_tibble %>%
	dplyr::mutate(contrast = myVar,
	Set = pathwayName)
})
names(gseas.PD) <- names(pathways)
do.call(rbind, gseas.PD) -> df
#Save results in table
df %>% filter(padj < pthres) %>% 
	arrange(desc(abs(NES)), padj) %>%
	print(n = Inf) %>%
	select(pathway, padj, NES) %>%
	mutate(padj =  round(padj, digit = 4),
	       NES = round(NES, digit = 3)) %>%
	write.table(., sep = ",", file = paste0("./", myVar,"_pathways.csv "), quote = F, row.names = F)
df %>% filter(padj < pthres) %>% 
	pull(leadingEdge) -> leadingEdge
names(leadingEdge) <- df %>% filter(padj < pthres) %>% pull(pathway)

#plot leadingEdge
unlist(leadingEdge) %>% table %>%
	sort %>%
	data.frame(gene_name = as.factor(names(.))) %>%
	arrange(Freq) %>%
	filter(Freq >= leadThres) %>%
	mutate(Type = case_when(
				grepl("^PSM", gene_name) ~ "Proteasome",
				grepl("^MR", gene_name) ~ "Mitochondrial ribosome",
				gene_name == "DAP3" ~ "Mitochondrial ribosome",	
				gene_name %in% c("COX7A2L", oxphos$gene_name) ~ "Mitochondrial resipratory chain",
				gene_name %in% c("AIFM1", "CS", "DLAT", "ECSIT", "IDH3A",
						 "IDH3B", "MDH1", "ME2", "FH", "PDHB",
						 "OGDH", 
						 "SUCLA2", "NNT") ~ "Mitochondrial enzyme/ associated",
				grepl("SLC", gene_name) ~ "Solute carrier",
				grepl("HSP", gene_name) ~ "Heatshock protein",
				grepl("SNC", gene_name) ~ "Synuclein",
				grepl("COL", gene_name) ~ "Collagens",
				grepl("KRT", gene_name) ~ "Keratine gene family",
				grepl("RPL", gene_name) ~ "Ribosomal",
				grepl("RPS", gene_name) ~ "Ribosomal",
				gene_name %in% (mgp %>%
						filter(cell_type == "Neurons") %>%
						pull(gene_name)) ~ "Neuron marker",
				grepl("^ATP6V", gene_name) ~ "V-type ATPases",

				gene_name %in% c("MARK1", "MAPK3", "MAP2K6") ~ "Kinase",
				gene_name %in% (mgp %>%
						filter(cell_type == "Oligodendrocytes") %>%
						pull(gene_name)) ~ "Oligo marker",
				gene_name %in% (mgp %>%
						filter(cell_type == "Astrocytes") %>%
						pull(gene_name)) ~ "Astrocyte marker",
				gene_name %in% (mgp %>%
						filter(cell_type == "Microglia") %>%
						pull(gene_name)) ~ "Microglia marker",
				gene_name %in% synapses ~ "Synapse",
				grepl("SYT", gene_name) ~ "Synapse",
				grepl("LAMA", gene_name) ~ "Laminin",
				grepl("LAMB", gene_name) ~ "Laminin",
				grepl("LAMP", gene_name) ~ "Lysosome a. mem. prot",
				grepl("^CCT", gene_name) ~ "Chaperonin containing CCT",
				TRUE ~ "Other")) %>%
	mutate(gene_name = factor(gene_name, levels = unique(gene_name))) -> df_leadingEdge
ggplot(df_leadingEdge, aes(x = gene_name, y = Freq, fill = Type)) +
	geom_bar(stat = "identity")+
#	coord_flip() +
	theme_cowplot() +
	labs(x = "Gene name", y = "Freq. of being in leading edge of enriched pathway") +
	theme(axis.text = element_text(size = 10)) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))	-> leplot

return(list(resTable = resM1, genePlot = gp, leadingEdge = leplot, resEnrich = df, leadingEdgeData = leadingEdge))
}

cellEstimates <- function(mat, df, celltype, highest = -0.15, lowest = -Inf, minus = F, mgp, rawdf, colBy = "Diagnosis"){
set.seed(123)
mgp %>% filter(cell_type == celltype) %>% pull(gene_name) -> cell
cell <- data.frame(gene_name = cell) %>% 
	dplyr::left_join(., 
			 rawdf %>% select(Accession, gene_name),
			 by = "gene_name",
			 multiple = "all") %>%
	filter(!is.na(Accession))
prcomp(log1p(mat[,colnames(mat) %in% cell$Accession]), scale = T, center = T)$rotation  %>%
	as.data.frame %>%
	mutate(Accession = rownames(.)) %>%
	filter(PC1 < highest, PC1 > lowest) %>% 
	pull(Accession) -> cell_II
cell_m <- log1p(mat[, colnames(mat) %in% cell_II])
colnames(cell_m) %>% data.frame(Accession = .) %>%
	dplyr::left_join(., rawdf, by = "Accession") %>%
	pull(gene_name) -> gn
colnames(cell_m) <- gn
prcomp(cell_m, scale = T, center = T) %>% 
	ggbiplot::ggbiplot(., group = df %>% pull(colBy)) +
  	scale_color_manual(values = cols_pd) +
	cowplot::theme_cowplot() -> p1
p1
prcomp(cell_m, scale = T, center = T) -> mypca
mypca$x %>%
	as.data.frame %>%
	select(PC1) %>%
	mutate(id = rownames(.))  -> estimate

#this is really bad...need to find a better conditioning
if(minus == FALSE) {
	estimate <- estimate %>%
		mutate(mgp = scales::rescale(PC1)) %>%
		select(id, mgp)
	} else {
		estimate <- estimate %>%
			mutate(mgp = 1-scales::rescale(PC1)) %>%
			select(id, mgp)
	}

df <- df %>% 
	dplyr::left_join(.,
			 estimate, by = "id")
df %>%
	ggplot(., aes(x = Diagnosis, y = mgp, col = Diagnosis)) +
	ggbeeswarm::geom_beeswarm() +
	geom_boxplot(fill = "transparent") +
	ggpubr::stat_compare_means() +
	ggpubr::stat_compare_means(comparisons = list(c("PSP", "Control"),
						      c("PSP", "PD"),
						      c("Control", "PD"),
						      c("MSA", "PD"),
						      c("MSA", "Control"),
						      c("MSA", "PSP"))) +
	scale_color_manual(values = cols_pd) +
	cowplot::theme_cowplot() +
	labs(y = paste0(celltype, " MGP"), tag ="A") -> p2

mypca$rotation[,1] %>% abs %>% sort -> topProteins
topProteins <- names(tail(topProteins, n = 3))
cell_m %>% as.data.frame %>%
	mutate(id = rownames(.)) %>%
	left_join(., df %>% 
			select(mgp, Diagnosis, id),
		by = "id") %>%
	reshape2::melt(., id.vars = c("id", "mgp", "Diagnosis")) %>%
	dplyr::rename(Accession = variable) %>%
	left_join(., rawdf %>% select(Accession, gene_name), by = "Accession") %>%
	filter(Accession %in% topProteins) %>%
	ggplot(., aes(x = Diagnosis, y = value, col = mgp)) +
	#ggbeeswarm::geom_beeswarm() +
	scale_color_viridis() +
	geom_jitter(height = 0, width = .1) +
	geom_boxplot(fill = "transparent", width = .4, outlier.shape = NA) +
	#geom_violin(fill = "transparent") +
	facet_wrap(~gene_name, nrow = 1) +
	labs(y = "Expression value", col = "MGP", tag = "C") +
	cowplot::theme_cowplot() -> p3
	p <- (p1 + p2)/ p3
	return(list("estimate" = estimate, "plot" = p, "loadings" = mypca$rotation)) 
}

runEnrich <- function(df, myGenes = NULL, thres = 0.5) {
library(parallel)
WebGestaltR::listGeneSet()[,1] -> dbs
dbs_selected <- dbs[c(1:11, 59:61,64)]
if (is.null(myGenes)) {
	myGenes <- df %>% filter(adj.P.Val < 0.05, abs(logFC) > thres) %>%
		    pull(gene_name)
}
print(myGenes)
ws <- mclapply(dbs_selected, function(db) {
print(paste0("ORA in: ", db))
up <- tryCatch({
	WebGestaltR(enrichDatabase = db,
     	    interestGene = myGenes, 
	    interestGeneType = "genesymbol",
	    referenceGene = df$gene_name,
	    referenceGeneType = "genesymbol",
	hostName="https://www.webgestalt.org/",
	isOutput = F,
	organism = "hsapiens")
	}, error = function(e) {
	NULL})		
 return(up)      
})
names(ws) <- dbs_selected
ws
 }


.toList <- function(data) {
	if (length(data)>2) {
		print(data)
#		data <- data[!is.na(data) && !is.null(data)]
		data <- data[!is.na(data) & !is.null(data)]
		# replace % in gene set names to _, because png treats % in filename specially
		data1 <- cbind(rep(gsub('%', '_', data[1], fixed=TRUE), length(data)-2), rep(data[2], length(data)-2), data[c(-1,-2)])
		return(data1)
	} else {
		return(NULL)
	}
}
readGmt <- function(gmtFile, cache=NULL) {
#####Change a gmt file to a three column matrix (gene set name, gene set description and genes)#######
	if (startsWith(gmtFile, "http://") || startsWith(gmtFile, "https://")) {
		response <- cacheUrl(gmtFile, cache)
		if (response$status_code == 200) {
			data <- unlist(strsplit(content(response, "text"), "\n", fixed=TRUE))
		} else {
			stop(webRequestError(response))
		}
	} else {
		if (file_ext(gmtFile) != "gmt") {
			stop(gmtFormatError("empty"))
		}
		# remove BOM in some windows files
		data <- gsub("\xEF\xBB\xBF", "", readLines(gmtFile, skipNul=TRUE), useBytes=TRUE)
	}
	data <- strsplit(data, "\t", useBytes=TRUE)
	data <- lapply(data,.toList)
	data <- do.call("rbind",data)
	if (is.null(data)) {
		stop(gmtFormatError("incorrect"))
	} else {
		data <- as.data.frame(data, stringsAsFactors=FALSE)
		colnames(data) <- c("geneSet", "description", "gene")
		return(data)
	}
}
readGMT <- readGmt

#' Load gene set data
#'
#' @inheritParams WebGestaltR
#'
#' @return A list of \code{geneSet}, \code{geneSetDes}, \code{geneSetDag}, \code{geneSetNet}, \code{standardId}.
#' \describe{
#'  \item{geneSet}{Gene set: A data frame with columns of "geneSet", "description", "genes"}
#'  \item{geneSetDes}{Description: A data frame with columns of two columns of gene set ID and description}
#'  \item{geneSetDag}{DAG: A edge list data frame of two columns of parent and child. Or a list of data frames if multilple databases are given.}
#'  \item{geneSetNet}{Network: A edge list data frame of two columns connecting nodes. Or a list of data frames if multilple databases are given.}
#'  \item{standardId}{The standard ID of the gene set}
#' }
#'
#' @importFrom dplyr select distinct filter %>%
#' @importFrom httr modify_url
#' @export
loadGeneSet <- function(organism="hsapiens", enrichDatabase=NULL, enrichDatabaseFile=NULL, enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL, cache=NULL, hostName="https://www.webgestalt.org/") {
	# TODO: multiple custom database ID types?
	geneSet <- NULL    ##gene sets
	geneSetDes <- NULL ##gene set description file
	geneSetDag <- list() ##gene set DAG file
	geneSetNet <- list() ##gene set network file
	standardId <- NULL

	if (organism != "others" && !is.null(enrichDatabaseFile) && is.null(enrichDatabaseType)) {
		stop("The ID type should be given in enrichDatabaseType for custom GMT files, e.g. genesymbol.")
	}
	# necessary because loop with length used below. enrichDatabase will skip NULL
	if (!is.vector(enrichDatabaseFile)) {
		enrichDatabaseFile = list(enrichDatabaseFile)
	}
	if (!is.vector(enrichDatabaseDescriptionFile)) {
		enrichDatabaseDescriptionFile = list(enrichDatabaseDescriptionFile)
	}
	if (length(enrichDatabaseFile) != length(enrichDatabaseDescriptionFile)) {
		stop("The number of custom database and its description files should be equal. Use NULL for placeholder.")
	}
	if (organism != "others") {  # supported organism
		geneSetInfo <- listGeneSet(organism=organism, hostName=hostName, cache=cache)
		#  load build-in databases
		for (enrichDb in enrichDatabase) {
			if (is.null(enrichDb) || enrichDb == "others") { next }  # just for backward compatibility
			if (!enrichDb %in% geneSetInfo$name) {
				warning("Database ", enrichDb, " is not supported")
				next
			}
			# get the ID type of the enriched database, such as entrezgene or phosphsiteSeq
			thisStandardId <- filter(geneSetInfo, .data$name==enrichDb)[[1, "idType"]]
			if (!is.null(standardId) && standardId != thisStandardId) {
				stop("Databases have inconsistent ID types. Mixed gene annotation databases with phosphosite databases?")
			}
			standardId <- thisStandardId

			#########Read GMT file from the existing database###########
			if (startsWith(hostName, "file://")) {
				gmtPath <- removeFileProtocol(file.path(hostName, "geneset", paste0(paste(organism, enrichDb, standardId, sep="_"), ".gmt")))
				thisGeneSet <- readGmt(gmtPath)
				thisGeneSet$database <- enrichDb # add a column for database source
				geneSet <- rbind(geneSet, thisGeneSet)
			} else {
				gmtUrl <- modify_url(file.path(hostName, "api", "geneset"), query=list(organism=organism, database=enrichDb, standardId=standardId, fileType="gmt"))
				thisGeneSet <-  readGmt(gmtUrl, cache=cache)
				thisGeneSet$database <- enrichDb
				oriCnt <- nrow(thisGeneSet)
				thisGeneSet <- thisGeneSet %>% filter(!(.data$geneSet %in% unique(!!geneSet$geneSet)))
				if (nrow(thisGeneSet) < oriCnt) {
				  warning(paste("Duplicate gene set names in", enrichDb, "have been ignored."))
				}
				geneSet <- rbind(geneSet, thisGeneSet)
			}

			#########Read the description file#############
			#geneSetDes <- rbind(geneSetDes, .loadGeneSetData(hostName, organism, enrichDb, standardId, "des", cache))
			thisGeneSetDes <- .loadGeneSetData(hostName, organism, enrichDb, standardId, "des", cache)
			if (!is.null(thisGeneSetDes) && !is.null(geneSetDes)) {
			  thisGeneSetDes <- thisGeneSetDes %>% filter(!(.data$geneSet %in% unique(!!geneSetDes$geneSet)))
			}
			geneSetDes <- rbind(geneSetDes, thisGeneSetDes)

			###########Try to load the DAG file#################
			# assignment considering possible return of NULL
			# list[] <- NULL will delete the element
			geneSetDag[enrichDb] <- list(.loadGeneSetData(hostName, organism, enrichDb, standardId, "dag", cache))

			###########Try to load the network file if the gene sets are generated from the network##########
			geneSetNet[enrichDb] <- list(.loadGeneSetData(hostName, organism, enrichDb, standardId, "net", cache))
		}

		# load local database files
		for (i in 1:length(enrichDatabaseFile)) {
			enrichDbFile <- enrichDatabaseFile[[i]]
			if (is.null(enrichDbFile)) { next }

			thisGeneSet <- idMapping(organism=organism, dataType="gmt", inputGeneFile=enrichDbFile, sourceIdType=enrichDatabaseType, targetIdType=NULL, mappingOutput=FALSE, cache=cache, hostName=hostName)
			thisStandardId <- thisGeneSet$standardId  # should be just enrichDatabaseType here
			if (!is.null(standardId) && standardId != thisStandardId) {
				stop("Databases have inconsistent ID types. Mixed gene annotation databases with phosphosite databases?")
			}
			standardId <- thisStandardId

			thisGeneSet <- thisGeneSet$mapped %>% select(.data$geneSet, .data$description, gene=standardId) %>% distinct()

			# load local description files
			enrichDbDesFile <- enrichDatabaseDescriptionFile[[i]]
			if (!is.null(enrichDbDesFile)) {
				thisGeneSetDes <- .loadEnrichDatabaseDescriptionFile(thisGeneSet, enrichDbDesFile)
				if (!is.null(thisGeneSetDes) && !is.null(geneSetDes)) {
				  thisGeneSetDes <- thisGeneSetDes %>% filter(!(.data$geneSet %in% unique(!!geneSetDes$geneSet)))
				}
				geneSetDes <- rbind(geneSetDes, thisGeneSetDes)
			}

			fileName <- gsub(".gmt", "", basename(enrichDbFile), fixed=TRUE)
			thisGeneSet$database <- fileName
			oriCnt <- nrow(thisGeneSet)
			thisGeneSet <- thisGeneSet %>% filter(!(.data$geneSet %in% unique(!!geneSet$geneSet)))
			if (nrow(thisGeneSet) < oriCnt) {
			  warning(paste("Duplicate gene set names in", fileName, "have been ignored."))
			}
			geneSet <- rbind(geneSet, thisGeneSet)

			geneSetDag[fileName] <- list(NULL) # correct way to assign NULL to list element
			geneSetNet[fileName] <- list(NULL)
		}
	} else { # custom organisms
		for (i in 1:length(enrichDatabaseFile)) {
			enrichDbFile <- enrichDatabaseFile[[i]]
			if (is.null(enrichDbFile)) { next }
			thisGeneSet <- readGmt(enrichDbFile)

			enrichDbDesFile <- enrichDatabaseDescriptionFile[[i]]
			if (!is.null(enrichDbDesFile)) {
				thisGeneSetDes <- .loadEnrichDatabaseDescriptionFile(thisGeneSet, enrichDbDesFile)
				if (!is.null(thisGeneSetDes) && !is.null(geneSetDes)) {
				  thisGeneSetDes <- thisGeneSetDes %>% filter(!(.data$geneSet %in% unique(!!geneSetDes$geneSet)))
				}
				geneSetDes <- rbind(geneSetDes, thisGeneSetDes)
			}
			fileName <- gsub(".gmt", "", basename(enrichDbFile), fixed=TRUE)
			thisGeneSet$database <- fileName
			oriCnt <- nrow(thisGeneSet)
			thisGeneSet <- thisGeneSet %>% filter(!(.data$geneSet %in% unique(!!geneSet$geneSet)))
			if (nrow(thisGeneSet) < oriCnt) {
				warning(paste("Duplicate gene set names in", enrichDb, "have been ignored."))
			}
			geneSet <- rbind(geneSet, thisGeneSet)

			geneSetDag[fileName] <- list(NULL)
			geneSetNet[fileName] <- list(NULL)
		}
	}
	if (is.null(geneSet)) { stop(enrichDatabaseError(type="empty")) }

	# unlist if just one database, for backward compatibility
	if (length(geneSetDag) == 1) { geneSetDag <- geneSetDag[[1]] }
	if (length(geneSetNet) == 1) { geneSetNet <- geneSetNet[[1]] }

	if (length(unique(geneSet$database)) == 1) {
		# remove database column for single source
		geneSet$database<- NULL
	}
	re <- list(geneSet=geneSet, geneSetDes=geneSetDes, geneSetDag=geneSetDag, geneSetNet=geneSetNet,standardId=standardId)
	return(re)
}

#' @importFrom httr content
#' @importFrom readr read_tsv
.loadGeneSetData <- function(hostName, organism, database, standardId, fileType, cache=NULL) {
	# read gene set files from API or returns NULL
	if (startsWith(hostName, "file://")) {
		geneSetPath <- removeFileProtocol(file.path(hostName, "geneset", paste(paste(organism, database, standardId, sep="_"), fileType, sep=".")))
		if (file.exists(geneSetPath)) {
			geneSetData <- read_tsv(geneSetPath, col_names=FALSE, col_types="cc", quote="")
		} else {
			geneSetData <- NULL
		}
	} else {
		geneSetUrl <- file.path(hostName,"api","geneset")
		response <- cacheUrl(geneSetUrl, cache=cache, query=list(organism=organism, database=database, standardId=standardId, fileType=fileType))
		if (response$status_code == 200) {
			geneSetData <- read_tsv(content(response), col_names=FALSE, col_types="cc", quote="")
		} else {
			geneSetData <- NULL
		}
	}
	if (!is.null(geneSetData) && fileType == "des") {
		colnames(geneSetData) <- c("geneSet", "description")
		geneSetData <- geneSetData %>% distinct(.data$geneSet, .keep_all=TRUE)
	}
	return(geneSetData)
}

#' @importFrom readr read_tsv
#' @importFrom tools file_ext
.loadEnrichDatabaseDescriptionFile <- function(geneSet, enrichDatabaseDescriptionFile) {
	if (file_ext(enrichDatabaseDescriptionFile) != "des") {
		warning(descriptionFileError("format"))
		return(NULL)
	} else {
		geneSetDes <- read_tsv(enrichDatabaseDescriptionFile, col_names=c("geneSet", "description"), col_types="cc", quote="") %>%
			distinct(geneSet, .keep_all=TRUE)
		if (ncol(geneSetDes)!=2) {
			warning(descriptionFileError("columnNum"))
			return(NULL)
		} else {
			if (length(intersect(unique(geneSet$geneSet), geneSetDes$geneSet)) < 0.6 * length(unique(geneSet[,1]))) {
				warning(descriptionFileError("overlap"))
				return(NULL)
			} else {
				return(geneSetDes)
			}
		}
	}
}
