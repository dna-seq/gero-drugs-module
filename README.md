### pgbprs

Code for gero-drugs module.

## Current description

_note that it is constantly changing_

Tool uses the *PharmGKB* (database)[https://www.pharmgkb.org/] as a source for information about influence of genetic variants on drugs by using of *polars* (package)[https://www.pola.rs/].


# File structure

- *inputdata* is a folder for the Source data: var\_drug\_ann.tsv table with pharmacological description of genetic varints (include critical - RSIDs), study\_parameters.tsv table with statistical parameters of effects (include critical - effect values) and (only for current state) toy-rsids.tsv table with test input (with RSIDs and Ref/Alt nucleotides) of variants.
- *tempdata* contains temporary table for genetic variant interpretation annotation\_tab.csv, this table is produced by merging and filtering of source tables from *inputdata* and contains only usable data. Subfolder *excluded* contains excluded datasets which are not suitable for current analysis due to different kinds of inconsistencies, flaws and incompleteness (see below).
- *output* is designed for an output report table, now there is a some small test output.
- PGBPRS.py3 is the main tool script, requires (for now) two tables from PharmGKB and some (not yet standardized) table input with RSIDs and corresponded nucleotides.
- PGBPRS.R is an investigational script for PharmGKB tables parsing and analysis, needed only for developmental reasons. According to it we decided to use var\_drug\_ann.tsv and study\_parameters.tsv tables.
- SharedCols.png is a technical output image for shared columns determination across all available PharmGKB tables.


## setting up

Install conda/micromamba environment
-------------------------
Annotation modules and dvc are included in the conda environment that goes with the project.
The environment can be setup either with Anaconda or micromamba (superior version of conda).
Here is installation instruction for micromamba:
```
wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
```
We can use ./micromamba shell init ... to initialize a shell (.bashrc) and a new root environment in ~/micromamba:
```
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
```
To create a micromamba environment use:
```
micromamba create -f environment.yaml
micromamba activate gero-drugs
```

The instructions above are provided for Linux and MacOS (note: in MacOS you have to install wget).
For Windows you can either install Linux Subsystem or use Windows version of anaconda.

Runing gero-drugs
--------------

To annotate samples you should use gero_drugs.py command line tool.
First, you should activate your environment and prepare the annotations by running
```
micromamba activate gero-drugs
python gero_drugs.py init
```
Then you can generate report by running:
```
python gero_drugs run inputdata/toy-rsids.tsv tempdata/annotation_tab.tsv
```
where you can substitute the default values with the ones you want.

# Workflow

_Initialization_ (actually these steps should be separated, now it isn't separated)

1. Read PharmGKB tables (in future we will use a normal SQL DB, not CSV reading)
2. Select only columns 'Variant Annotation ID', 'Variant/Haplotypes', 'Drug(s)', 'Phenotype Category', 'Significance', 'Sentence' for var\_drug\_ann.tsv table and only 'Variant Annotation ID', 'Allele Of Frequency In Cases', 'Allele Of Frequency In Controls', 'P Value', 'Ratio Stat Type', 'Ratio Stat', 'Confidence Interval Start', 'Confidence Interval Stop' columns for study\_parameters.tsv (fields may change during development)
3. Join both tables to the raw annotationTab table by share field 'Variant Annotation ID'.
4. Exclude inappropriate data entries (which cannot be processed now without extra preparation) and save to separate 'Problematic' tables, detailed description is below.
5. Drop entries without available effect value ('Ratio Stat').
6. Save the ready annotationTab table.

_Interpretation_ (remind that this section is too raw)

1. Read input with genotype (currently it is a simplistic table with RSID and nucleotides, in future it should be some sort of VCF).
2. Join input and annotationTab tables by RSID.
3. Save produced report (currently future processing is not yet developed, this table may be transformed to some pivot table in an interactive viewer).


# Data inconsistencies and Problematic tables

The PharmGKB database isn't just database actually, it is knowledge base and significant part is not fully formalized - not transformed to ordered numbers and parameters by strict rules. By this reason a lot of entries require additional transformations, development of new conversion rules and manual curation. These datasets called Problematic tables and its will be added in future after development of new conversion rules. There is the current list of excluded data categories and open questions:

- The var\_drug\_ann in it's main ID column 'Variant/Haplotypes' sometimes contains multiple entries (written by comma) in one line, these cases have no known ways for handling and saved to problematicMultiVarHapEntries.csv table.
- The var\_drug\_ann has similar problem with 'Drug(s)' - multiple entries in one line, these lines saved to problematicMultiDrugEntries.csv table.
- The same problem has for 'Phenotype Category' column (instead of only three possible values) with multiple values, saved to problematicMultiPhenoCat.csv table.
- Sometimes Drug isn't an exact substance but drug class, this is a challenge for interpretation, all these cases are excluded (by simple language trick) and saved to problematicDrugClass.csv table.
- Many entries in PharmGKB not about SNPs (as RSIDs) but about haplotypes (Cytochromes, HLA and etc.) and it't need additional conversion or different way of analysis (we have an idea with additional correspondence table), saved to problematicHaps.csv table.
- Sometimes the effect value in study\_parameters.tsv is presented but 'Ratio Stat Type', statistical category of value (RR, OR or HR) is unknown literally, these cases are excluded to separate table problematicUnknownStatType.csv for manual curation.
- For interpretation we need for each SNP to know what is Control nucleotide and Sample nucleotide in each study (there is a problem that Alt and Red may alternate during progess in genomics and an old reference may became an alternative) simultaneously with information about actual reference and "patient's" genotype, i.e. we need to know Ref and Alt from each study. Many entries have no available Ref and/or Alt, we excluded and saved these entries to problematicUnclearRefAltNucl.csv tables for manual curation.
- The PRS may be constructed from numerical relations of effects and variants but many entries in study\_parameters.tsv have no value (have only null) even instead of p-value presence. We can't use null for numerical analysis and all these cases just dropped.
- By all these reasons we have only small ~140-line pharmacogemonic table.

# Unresolved questions

- A lot of data can't be used due to additional processing / manual curation requirements.
- Table input/output instead of SQL and Web.
- Not a OpenCravat/OakVar module yet, just a bunch of Python code.
- It is unknown how to add support of heterozygosity.
- Actual interpretation is not added yet even substitution of nucleotides from input.
 

