
# get file name list
files <- list.files(path = 'pgbprs', pattern = '.tsv')
# init list of dataframes
tables <- list()
# go over files
for (name in files) {
	print(name)
	# read table
	temp <- read.table(paste0('pgbprs/', name), header = T, sep = '\t', quote = "", fill = T, stringsAsFactors = F)
	# convert colnames to lowercase
	colnames(temp) <- tolower(colnames(temp))
	# remove final ".s." from gene.s and phenotype.s
	colnames(temp) <- gsub(pattern = ".s.$", replacement = "", x = colnames(temp))
	# assign name for list
	tempname <- gsub('.tsv', '', name)
	# save dataframe to list
	tables[[tempname]] <- temp
}

# init dataframe for overlaps
sharedata <- data.frame(table1 = character(), table2 = character(), nt1 = integer(), nt2 = integer(), count = integer(), share = character(), stringsAsFactors = F)
# count share columns across tables
# note that this way is blind for table "Relationships" because it contains Chemical, Disease, Gene, Haplotype and Variant in columns entity1_type and entity2_type
for (i in 1:length(tables)) {
	for (j in 1:length(tables)) {
		# get i and j roww
		irow <- colnames(tables[[i]])
		jrow <- colnames(tables[[j]])
		# get names
		iname <- names(tables)[i]
		jname <- names(tables)[j]
		# count overlap
		countval <- sum(irow %in% jrow)
		# get share columns
		shared <- irow[irow %in% jrow]
		# turn zero length to ''
		if (length(shared) == 0) { shared <- ''}
		# stitch multiple columns
		if (length(shared) > 1) { shared <- paste0(shared, collapse = ',') }
		
		sharedata[nrow(sharedata) + 1, ] <- list(iname, jname, i, j, countval, shared)
	}
}

# plot shared "table"
png(filename = 'pgbprs/SharedCols.png', width = 600, height = 600, res = 100)
par(mar = c(8, 8, 0, 0))
plot(0, 0, type = 'n', axes = F, xlim = c(0, 16), ylim = c(0, 16), xlab = '', ylab = '')
text(x = sharedata$nt1, y = sharedata$nt2, labels = sharedata$count)
abline(a = 0, b = 1)
mtext(text = names(tables), side = 1, at = 1:15, las = 2)
mtext(text = names(tables), side = 2, at = 1:15, las = 2)
dev.off()

## Data extraction
# We start from "study_parameters" table (number 11) - only this table contains effect values for genetic variants
subsharedata <- sharedata[sharedata$table1 == 'study_parameters', ]
subsharedata <- subsharedata[subsharedata$table2 != 'study_parameters', ]
subsharedata <- subsharedata[subsharedata$count > 0, ]
# this table has only one shared column: 'variant.annotation.id' but across 3 tables - "var_drug_ann", "var_fa_ann' and "var_pheno_ann"
# var_pheno_ann.tsv: Contains associations in which the variant affects a phenotype, with or without drug information.
# var_drug_ann.tsv: Contains associations in which the variant affects a drug dose, response, metabolism, etc
# var_fa_ann.tsv: Contains in vitro and functional analysis-type associations.
# According to this - var_fa_ann should be excluded because it cames from non clinical data, var_pheno_ann contains ambiguous data (drug info may be or may be not presented), let's exclude this table. Focus on var_drug_ann only:
# 'study_parameters' <-- variant.annotation.id --> 'var_drug_ann'
# and "var_drug_ann" table has 'variant.haplotypes' column which is RSID i.e. this is an Input table for RSID however frequently 'variant.haplotypes' column has haplotype (usually CYP) instead of RSID or has multiple entries for one 'varian.annotation.id'.

## estimations
# number of single RSID entries
length(tables[['var_drug_ann']]$variant.haplotypes) # 11109 totally
temp <- grep(x = tables[['var_drug_ann']]$variant.haplotypes, pattern = '^rs', value = T) # RSID only - gsub by RE '^rs'
temp <- grep(x = temp, pattern = ',', value = T, invert = T) # single RSID only - !gsub by RE ','
length(temp) # 8592 totally but only 
length(unique(temp)) # ... 3415 only unique
# mean number of study parameters entry per RSID
mat <- matrix(0, nrow = length(unique(temp)), ncol = 2)
for (k in 1:length(unique(temp))) {
	rsid <- unique(temp)[k]
	# go over 'var_drug_ann'
	subtab <- tables[['var_drug_ann']][tables[['var_drug_ann']]$variant.haplotypes == rsid, ]
	# get variant annotation ID
	varannid <- subtab$variant.annotation.id
	# query variant annotation ID to 'study_parameters'
	subdat <- tables[['study_parameters']][tables[['study_parameters']]$variant.annotation.id %in% varannid, ]
	mat[k, 1] <- sum(temp == rsid)
	mat[k, 2] <- nrow(subdat)
}


