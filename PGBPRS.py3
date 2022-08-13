import polars as pl

# loading of tables
var_drug_ann = pl.read_csv("pgbprs/inputdata/var_drug_ann.tsv", has_header = True, sep = "\t")
study_parameters = pl.read_csv("pgbprs/inputdata/study_parameters.tsv", has_header = True, sep = "\t)

### RSID as an example
rsid = 'rs3745274'

# get var_drug_ann data by RSID, main value is 'Variant Annotation ID' as ID (VAID) and select only needed columns ('Variant/Haplotypes', 'Variant Annotation ID', 'Drug(s)', 'Phenotype Category', 'Significance')
annotationTab = var_drug_ann.filter(pl.col('Variant/Haplotypes') == rsid).select(['Variant/Haplotypes', 'Variant Annotation ID', 'Drug(s)', 'Phenotype Category', 'Significance'])

# make shrinked study_parameters DF
study_parameters_shrink = study_parameters.select(['Variant Annotation ID', 'P Value', 'Ratio Stat Type', 'Ratio Stat', 'Confidence Interval Start', 'Confidence Interval Stop'])

# merge variant data annotationTab with shrinked study_parameters table - may be used for SNP effect accounting
temp = annotationTab.join(study_parameters_shrink, on = 'Variant Annotation ID')


## Estimations
# merge shrinked study_parameters with var_drug_ann table
temp2 = var_drug_ann.join(study_parameters_shrink, on = 'Variant Annotation ID')
# remove haplotypes
temp2 = temp2[temp2['Variant/Haplotypes'].str.contains("^rs")]
# ~8300 kept
temp2 = temp2[temp2['Ratio Stat'].is_not_null()]
# ~1500 entries have non-null Ratio Stat value
temp2['Variant/Haplotypes'].value_counts()
# and ~880 have unique SNP
temp2['Drug(s)'].value_counts()
# ~200 unique drgus but estimation is very rough due to frequent multiple drug entries


