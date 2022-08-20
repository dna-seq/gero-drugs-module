import polars as pl

### Initial preparation of data - it is needed only once

## loading of source tables
# var_drug_ann is Variant-Drug table
var_drug_ann = pl.read_csv("pgbprs/inputdata/var_drug_ann.tsv", has_header = True, sep = "\t")
# study_parameters is Study Parameters table
study_parameters = pl.read_csv("pgbprs/inputdata/study_parameters.tsv", has_header = True, sep = "\t")

# ### RSID as an example
# rsid = 'rs3745274'

# select only needed columns from Variant-Drug table: 'Variant/Haplotypes', 'Variant Annotation ID', 'Drug(s)', 'Phenotype Category', 'Significance' and 'Sentence' for future annotation tab
var_drug_ann_shrink = var_drug_ann.select(['Variant Annotation ID', 'Variant/Haplotypes', 'Drug(s)', 'Phenotype Category', 'Significance', 'Sentence'])

# select only needed columns from Study Parameters table: 'Variant Annotation ID', 'P Value', 'Ratio Stat Type', 'Ratio Stat', 'Confidence Interval Start', 'Confidence Interval Stop' 
study_parameters_shrink = study_parameters.select(['Variant Annotation ID', 'Allele Of Frequency In Cases', 'Allele Of Frequency In Controls', 'P Value', 'Ratio Stat Type', 'Ratio Stat', 'Confidence Interval Start', 'Confidence Interval Stop'])

# merge shrinked Variant-Drug table with shrinked study_parameters table for assembly of Annotation table
annotationTab = var_drug_ann_shrink.join(study_parameters_shrink, on = 'Variant Annotation ID')

## Filtration
# ! exclude and save entries with multiple 'Variant/Haplotypes' values (contains ',')
problematicMultiVarHapEntries = annotationTab.filter(pl.col('Variant/Haplotypes').str.contains(","))
annotationTab = annotationTab.filter(pl.col('Variant/Haplotypes').str.contains(",") == False)
# ! exclude and save entries with multiple 'Drug(s)' values (contains ',' or '/') !
problematicMultiDrugEntries = annotationTab.filter(pl.col('Drug(s)').str.contains(",|/"))
annotationTab = annotationTab.filter(pl.col('Drug(s)').str.contains(",|/") == False)
# ! exclude and save entries with multiple 'Phenotype Category'
problematicMultiPhenoCat = annotationTab.filter(pl.col('Phenotype Category').str.contains(",|/"))
annotationTab = annotationTab.filter(pl.col('Phenotype Category').str.contains(",|/") == False)
# ! exclude and save entries with drug class instead of ecact class in 'Drug(s)' values (contains 's$')
# ...except 'us$' because tacrolimus and sirolimus are not classes
problematicDrugClass = annotationTab.filter(pl.col('Drug(s)').str.contains("[^u]s$"))
annotationTab = annotationTab.filter(pl.col('Drug(s)').str.contains("[^u]s$") == False)
# ! exclude and save entries with Haplotypes in 'Variant/Haplotypes'
problematicHaps = annotationTab.filter(pl.col('Variant/Haplotypes').str.contains("^rs") == False)
annotationTab = annotationTab.filter(pl.col('Variant/Haplotypes').str.contains("^rs"))
# ! exclude and save entries with "unknown" 'Ratio Stat Type'
problematicUnknownStatType = annotationTab.filter(pl.col('Ratio Stat Type').str.contains("Unknown"))
annotationTab = annotationTab.filter(pl.col('Ratio Stat Type').str.contains("Unknown") == False)
# ! exclude and save entries with null Ref and/or Alt nucleotudes in 'Allele Of Frequency In Cases' and 'Allele Of Frequency In Controls'
problematicUnclearRefAltNucl = annotationTab.filter((pl.col('Allele Of Frequency In Cases') == '') | (pl.col('Allele Of Frequency In Controls') == ''))
annotationTab = annotationTab.filter((pl.col('Allele Of Frequency In Cases') != '') & (pl.col('Allele Of Frequency In Controls') != ''))

# drop entries with null Ratio Stai - with unknown effect value and useless for PRS construction
annotationTab = annotationTab.filter(pl.col('Ratio Stat').is_not_null())

# save data with problematic entries
problematicMultiVarHapEntries.to_csv('pgbprs/tempdata/excluded/problematicMultiVarHapEntries.csv', sep = "\t")
problematicMultiDrugEntries.to_csv('pgbprs/tempdata/excluded/problematicMultiDrugEntries.csv', sep = "\t")
problematicMultiPhenoCat.to_csv('pgbprs/tempdata/excluded/problematicMultiPhenoCat.csv', sep = "\t")
problematicDrugClass.to_csv('pgbprs/tempdata/excluded/problematicDrugClass.csv', sep = "\t")
problematicHaps.to_csv('pgbprs/tempdata/excluded/problematicHaps.csv', sep = "\t")
problematicUnknownStatType.to_csv('pgbprs/tempdata/excluded/problematicUnknownStatType.csv', sep = "\t")
problematicUnclearRefAltNucl.to_csv('pgbprs/tempdata/excluded/problematicUnclearRefAltNucl.csv', sep = "\t")

## save filtered data to CSV
annotationTab.to_csv('pgbprs/tempdata/annotation_tab.csv', sep = "\t")


## Genotype analysis

# # load Anton's VCF as sample
# personvcf = pl.read_csv("antonkulaga.vcf", has_header = False, sep = "\t", comment_char = "#")

# load toy sample
variantsTab = pl.read_csv("pgbprs/inputdata/toy-rsids.tsv", has_header = True, sep = "\t", comment_char = "#")

# report table
reportTab = annotationTab.join(variantsTab, on = 'Variant/Haplotypes')

# save report table
reportTab.to_csv('pgbprs/output/report.csv', sep = "\t")


