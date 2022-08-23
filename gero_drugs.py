from beartype import beartype
import polars as pl
import click
from pathlib import Path
from typing import Union

@beartype
def prepare_annotations(base: Path):

    ### Initial preparation of data - it is needed only once
    inputdata = (base / "inputdata").absolute().resolve()
    print(f"initial data preparation started, input folder is {str(inputdata)}")

    ## loading of source tables
    # var_drug_ann is Variant-Drug table
    var_drug_ann = pl.read_csv(f"{str(inputdata)}/var_drug_ann.tsv", has_header = True, sep = "\t")
    # study_parameters is Study Parameters table
    study_parameters = pl.read_csv(f"{str(inputdata)}/study_parameters.tsv", has_header = True, sep = "\t")
    
    # ### RSID as an example
    # rsid = "rs3745274"
    
    # select only needed columns from Variant-Drug table: "Variant/Haplotypes", "Variant Annotation ID", "Drug(s)", "Phenotype Category", "Significance" and "Sentence" for future annotation tab
    var_drug_ann_shrink = var_drug_ann.select(["Variant Annotation ID", "Variant/Haplotypes", "Drug(s)", "Phenotype Category", "Significance", "Sentence"])
    
    # select only needed columns from Study Parameters table: "Variant Annotation ID", "P Value", "Ratio Stat Type", "Ratio Stat", "Confidence Interval Start", "Confidence Interval Stop" 
    study_parameters_shrink = study_parameters.select(["Variant Annotation ID", "Allele Of Frequency In Cases", "Allele Of Frequency In Controls", "P Value", "Ratio Stat Type", "Ratio Stat", "Confidence Interval Start", "Confidence Interval Stop"])
    
    # merge shrinked Variant-Drug table with shrinked study_parameters table for assembly of Annotation table
    annotation_tab = var_drug_ann_shrink.join(study_parameters_shrink, on="Variant Annotation ID")

    return filter_annotations(annotation_tab, base)


@beartype
def filter_annotations(annotation_tab: pl.DataFrame, base: Path) -> pl.DataFrame:
    assert base.exists(), "base path should exist!"
    ## Filtration
    # ! exclude and save entries with multiple "Variant/Haplotypes" values (contains ",")
    problematic_multi_var_hap_entries = annotation_tab.filter(pl.col("Variant/Haplotypes").str.contains(","))
    annotation_tab = annotation_tab.filter(pl.col("Variant/Haplotypes").str.contains(",") == False)
    # ! exclude and save entries with multiple "Drug(s)" values (contains "," or "/") !
    problematic_multi_drug_entries = annotation_tab.filter(pl.col("Drug(s)").str.contains(",|/"))
    annotation_tab = annotation_tab.filter(pl.col("Drug(s)").str.contains(",|/") == False)
    # ! exclude and save entries with multiple "Phenotype Category"
    problematic_multi_pheno_cat = annotation_tab.filter(pl.col("Phenotype Category").str.contains(",|/"))
    annotation_tab = annotation_tab.filter(pl.col("Phenotype Category").str.contains(",|/") == False)
    # ! exclude and save entries with drug class instead of ecact class in "Drug(s)" values (contains "s$")
    # ...except "us$" because tacrolimus and sirolimus are not classes
    problematic_drug_class = annotation_tab.filter(pl.col("Drug(s)").str.contains("[^u]s$"))
    annotation_tab = annotation_tab.filter(pl.col("Drug(s)").str.contains("[^u]s$") == False)
    # ! exclude and save entries with Haplotypes in "Variant/Haplotypes"
    problematic_haps = annotation_tab.filter(pl.col("Variant/Haplotypes").str.contains("^rs") == False)
    annotation_tab = annotation_tab.filter(pl.col("Variant/Haplotypes").str.contains("^rs"))
    # ! exclude and save entries with "unknown" "Ratio Stat Type"
    problematic_unknown_stat_type = annotation_tab.filter(pl.col("Ratio Stat Type").str.contains("Unknown"))
    annotation_tab = annotation_tab.filter(pl.col("Ratio Stat Type").str.contains("Unknown") == False)
    # ! exclude and save entries with null Ref and/or Alt nucleotudes in "Allele Of Frequency In Cases" and "Allele Of Frequency In Controls"
    problematic_unclear_ref_alt_nucl = annotation_tab.filter(
        (pl.col("Allele Of Frequency In Cases").is_null()) | (pl.col("Allele Of Frequency In Controls").is_null()))
    annotation_tab = annotation_tab.filter(
        (pl.col("Allele Of Frequency In Cases").is_not_null()) & (pl.col("Allele Of Frequency In Controls").is_not_null()))
    # ! exclude and save entries with the same Ref and Alt nucleotudes in "Allele Of Frequency In Cases" and "Allele Of Frequency In Controls"
    problematic_same_ref_alt_nucl = annotation_tab.filter(pl.col("Allele Of Frequency In Cases") == pl.col("Allele Of Frequency In Controls"))
    annotation_tab = annotation_tab.filter(pl.col("Allele Of Frequency In Cases") != pl.col("Allele Of Frequency In Controls"))

    # save data with problematic entries
    excluded = base / "tempdata" / "excluded"
    print(f"writing excluded data to {str(excluded)}")
    excluded.mkdir(exist_ok=True)
    problematic_multi_var_hap_entries.write_csv(f"{str(excluded)}/problematic_multi_var_hap_entries.tsv", sep = "\t")
    problematic_multi_drug_entries.write_csv(f"{str(excluded)}/problematic_multi_drug_entries.tsv", sep = "\t")
    problematic_multi_pheno_cat.write_csv(f"{str(excluded)}/problematic_multi_pheno_cat.tsv", sep = "\t")
    problematic_drug_class.write_csv(f"{str(excluded)}/problematic_drug_class.tsv", sep = "\t")
    problematic_haps.write_csv(f"{str(excluded)}/problematic_haps.csv", sep = "\t")
    problematic_unknown_stat_type.write_csv(f"{str(excluded)}/problematic_unknown_stat_type.tsv", sep = "\t")
    problematic_unclear_ref_alt_nucl.write_csv(f"{str(excluded)}/problematic_unclear_ref_alt_nucl.tsv", sep = "\t")
    problematic_same_ref_alt_nucl.write_csv(f"{str(excluded)}/problematic_same_ref_alt_nucl.tsv", sep = "\t")
    # drop entries with null Ratio Stat - with unknown effect value and useless for PRS construction
    annotation_tab_filtered = annotation_tab.filter(pl.col("Ratio Stat").is_not_null())

    filter_annotations_path = base / "tempdata" / "annotation_tab.tsv"
    annotation_tab_filtered.write_csv(f"{str(filter_annotations_path)}", sep = "\t")
    print(f"filtered annotations saved to {filter_annotations_path}")
    return annotation_tab_filtered

@beartype
def analyze(sample: str, annotation_tab, base: Path) -> Union[str, Path]:
    ## Genotype analysis
    # # load Anton"s VCF as sample
    # personvcf = pl.read_csv("antonkulaga.vcf", has_header = False, sep = "\t", comment_char = "#")
    # load toy sample
    variants_tab = pl.read_csv(sample, has_header=True, sep="\t", comment_char="#")
    # prepare report table
    report_tab = annotation_tab.join(variants_tab, on="Variant/Haplotypes")
    # create a separate column for output
    report_tab = report_tab.with_column(pl.col('Ratio Stat').alias('Effect'))
    # perform interprepation - substitute 'Alt' from genome to 'Allele Of Frequency In Cases' to table, it they don't match - invert effect vale
    report_tab = report_tab.with_column(pl.when(pl.col('Allele Of Frequency In Cases') == pl.col('alt')).then(pl.col('Effect')).otherwise((1/pl.col('Effect')).round(3)))
    # remove useless now columns
    report_tab = report_tab.drop(['Variant Annotation ID', 'Ratio Stat', 'Confidence Interval Start', 'Confidence Interval Stop', 'P Value', 'alt', 'ref'])
    # save report table
    report_path = f"{str(base)}/output/report.tsv"
    report_tab.write_csv(report_path, sep="\t")
    print(f"successfully wrote report to {report_path}")
    return report_path

@click.group()
def app():
    print("running gero-drugs application")

@app.command("init")
@click.argument("base", default=".")
def init(base: str = "."):
    annotation_tab = prepare_annotations(Path(base).absolute().resolve())
    return annotation_tab

@app.command("run")
@click.argument("sample", type=click.Path(exists=True), default="./inputdata/toy-rsids.tsv")
@click.argument("annotations", type=click.Path(exists=True), default="./tempdata/annotation_tab.tsv")
@click.argument("base", default=".")
def run(sample: str, annotations: str, base: str):
    base_folder = Path(base).absolute().resolve()
    print(f"initializing the annotations with basefolder ${str(base_folder)}")
    annotation_tab = pl.read_csv(annotations, sep="\t", has_header=True, comment_char="#")
    return analyze(sample, annotation_tab, base_folder)


if __name__ == "__main__":
    app()