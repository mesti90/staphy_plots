from pathlib import Path
import pandas as pd

configfile:"config/config.yaml"

OUT=Path(config["outdir"])
TABLEDIR=OUT/"{interval}"/"tables"
FILTERDIR=OUT/"{interval}"/"filtered"

DATASETS={
	"all":config["input"]["all"],
	"mrsa":config["input"]["mrsa"],
}

DATASET_NAMES=list(DATASETS.keys())

REGIONS=sorted(config["region_to_continent"].keys())
CONTINENTS=sorted(set(config["region_to_continent"].values()))

REGION_FILE_NAMES={region:region.replace(" ","_") for region in REGIONS}
REGION_FILE_TO_NAME={value:key for key,value in REGION_FILE_NAMES.items()}

ST_COLUMN=config["columns"]["mlst"]
_st_df=pd.read_csv(DATASETS["all"],sep="\t",usecols=[ST_COLUMN],dtype=str)
_st_counts=_st_df[ST_COLUMN].dropna().value_counts()

TAR_COLORS={
	"none":"#808080",
	"tarM":"#e6194b",
	"tarM+tarP":"#3cb44b",
	"tarM+tarP+tarS":"#ffe119",
	"tarM+tarS":"#4363d8",
	"tarP":"#f58231",
	"tarP+tarS":"#911eb4",
	"tarS":"#46f0f0",
}

ST_PALETTE=[
	"#e6194b","#3cb44b","#ffe119","#4363d8","#f58231",
	"#911eb4","#46f0f0","#f032e6","#bcf60c","#fabebe",
	"#008080","#e6beff","#9a6324","#800000","#aaffc3",
	"#808000","#ffd8b1","#000075",
	"#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
	"#ffff33","#a65628","#f781bf","#1b9e77","#d95f02",
	"#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d",
	"#1f78b4","#33a02c","#fb9a99","#6a3d9a",
]

ST_COLOR_MAP = {}
for i, st in enumerate(_st_counts.index):
	mycolor = ST_PALETTE[i%len(ST_PALETTE)]
	ST_COLOR_MAP[st] = mycolor
	st = float(st)
	ST_COLOR_MAP[st] = mycolor
	st = int(st)
	ST_COLOR_MAP[st] = mycolor
	st = str(st)
	ST_COLOR_MAP[st] = mycolor
	ST_COLOR_MAP[f"ST{int(st)}"] = mycolor
	
MLST_TAR_COLORS=[
	"#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
	"#ffff33","#a65628","#f781bf","#999999","#66c2a5",
	"#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f",
	"#e5c494","#1b9e77","#d95f02","#7570b3","#e7298a",
	"#66a61e","#e6ab02","#a6761d","#1f78b4","#33a02c",
	"#fb9a99","#6a3d9a","#b2df8a","#cab2d6","#ffed6f",
]

REGION_COLORS=[
	"#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
	"#ffff33","#a65628","#f781bf","#1b9e77","#d95f02",
	"#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d",
	"#1f78b4","#33a02c","#fb9a99","#6a3d9a","#b15928",
]

CONTINENT_COLORS=[
	"#e6194b","#3cb44b","#ffe119","#4363d8","#f58231",
	"#911eb4","#46f0f0","#f032e6","#bcf60c","#fabebe",
	"#008080","#e6beff","#9a6324","#800000","#aaffc3",
	"#808000","#ffd8b1","#000075",
]

REGION_COLOR_MAP={
	region:REGION_COLORS[i%len(REGION_COLORS)]
	for i,region in enumerate(REGIONS)
}

CONTINENT_COLOR_MAP={
	continent:CONTINENT_COLORS[i%len(CONTINENT_COLORS)]
	for i,continent in enumerate(CONTINENTS)
}

WORLD_COLOR_MAP={"world":"#377eb8"}

YEAR_INTERVALS=[
	f"{start}-{end}"
	for start,end in config["year_intervals"]
]

def get_interval(wc):
	for start,end in config["year_intervals"]:
		if wc.interval==f"{start}-{end}":
			return start,end
	raise ValueError(f"Unknown interval: {wc.interval}")


rule all:
	input:
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_heatmap.{dataset}.{interval}.region.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_heatmap.{dataset}.{interval}.continent.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"serotop.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"serotop.{dataset}.{interval}.region.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"serotop.{dataset}.{interval}.continent.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_rarefaction.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_rarefaction.{dataset}.{interval}.region.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_rarefaction.{dataset}.{interval}.continent.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_prevalence_pie.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_prevalence_pie.{dataset}.{interval}.region.{region_file}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			region_file=REGION_FILE_NAMES.values(),
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_prevalence_pie.{dataset}.{interval}.continent.{continent}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			continent=CONTINENTS,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_over_time.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_over_time.{dataset}.{interval}.region.{region_file}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			region_file=REGION_FILE_NAMES.values(),
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_over_time.{dataset}.{interval}.continent.{continent}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			continent=CONTINENTS,
		),
		expand(
			OUT/"{interval}"/"all"/"mrsa_over_time.all.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
		),
		expand(
			OUT/"{interval}"/"all"/"mrsa_over_time.all.{interval}.region.{region_file}.pdf",
			interval=YEAR_INTERVALS,
			region_file=REGION_FILE_NAMES.values(),
		),
		expand(
			OUT/"{interval}"/"all"/"mrsa_over_time.all.{interval}.continent.{continent}.pdf",
			interval=YEAR_INTERVALS,
			continent=CONTINENTS,
		),
		OUT/"stats/serotype_logistic_regression_all.xlsx",
		expand(
			OUT/"{interval}"/"{dataset}"/"serotype_logistic_regression.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"tar_over_time.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"tar_over_time.{dataset}.{interval}.region.{region_file}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			region_file=REGION_FILE_NAMES.values(),
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"tar_over_time.{dataset}.{interval}.continent.{continent}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			continent=CONTINENTS,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_tar_prevalence_pie.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_tar_prevalence_pie.{dataset}.{interval}.region.{region_file}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			region_file=REGION_FILE_NAMES.values(),
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_tar_prevalence_pie.{dataset}.{interval}.continent.{continent}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			continent=CONTINENTS,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_tar_over_time.{dataset}.{interval}.world.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_tar_over_time.{dataset}.{interval}.region.{region_file}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			region_file=REGION_FILE_NAMES.values(),
		),
		expand(
			OUT/"{interval}"/"{dataset}"/"mlst_tar_over_time.{dataset}.{interval}.continent.{continent}.pdf",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
			continent=CONTINENTS,
		),


rule prepare_dataset:
	input:
		tsv=lambda wc:DATASETS[wc.dataset]
	output:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	params:
		date_column=config["columns"]["date"],
		region_column=config["columns"]["region"],
		continent_column=config["columns"]["continent"],
		minyear=lambda wc:get_interval(wc)[0],
		maxyear=lambda wc:get_interval(wc)[1],
		region_to_continent=config["region_to_continent"]
	script:
		"scripts/prepare_dataset.py"


rule mlst_summary:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=OUT/"{interval}"/"{dataset}_mlst_summary.tsv"
	params:
		mlst_column=config["columns"]["mlst"],
		region_column=config["columns"]["region"],
		continent_column=config["columns"]["continent"],
		mrsa_column=config["mrsa"]["column"],
		mrsa_value=config["mrsa"]["value"]
	script:
		"scripts/mlst_summary.py"


### RAREFACTION ###

rule mlst_rarefaction_world:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_mlst_rarefaction_world.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_rarefaction.{dataset}.{interval}.world.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="world",
		geography_column=None,
		step=config["rarefaction"]["step"],
		repeats=config["rarefaction"]["repeats"],
		seed=config["rarefaction"]["seed"],
		geography_color_map=WORLD_COLOR_MAP
	script:
		"scripts/mlst_rarefaction.py"


rule mlst_rarefaction_region:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_mlst_rarefaction_region.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_rarefaction.{dataset}.{interval}.region.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_color_map=REGION_COLOR_MAP,
		step=config["rarefaction"]["step"],
		repeats=config["rarefaction"]["repeats"],
		seed=config["rarefaction"]["seed"]
	script:
		"scripts/mlst_rarefaction.py"


rule mlst_rarefaction_continent:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_mlst_rarefaction_continent.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_rarefaction.{dataset}.{interval}.continent.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		step=config["rarefaction"]["step"],
		repeats=config["rarefaction"]["repeats"],
		seed=config["rarefaction"]["seed"],
		geography_color_map=CONTINENT_COLOR_MAP
	script:
		"scripts/mlst_rarefaction.py"


### MLST OVER TIME ###

rule mlst_over_time_world:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_over_time_world.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_over_time.{dataset}.{interval}.world.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="world",
		geography_column=None,
		geography_value=None,
		threshold=config["timeseries"]["threshold"],
		st_color_map=ST_COLOR_MAP
	script:
		"scripts/mlst_over_time.py"


rule mlst_over_time_region:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_over_time_region_{region_file}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_over_time.{dataset}.{interval}.region.{region_file}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_value=lambda wc:REGION_FILE_TO_NAME[wc.region_file],
		threshold=config["timeseries"]["threshold"],
		st_color_map=ST_COLOR_MAP
	script:
		"scripts/mlst_over_time.py"


rule mlst_over_time_continent:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_over_time_continent_{continent}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_over_time.{dataset}.{interval}.continent.{continent}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_value=lambda wc:wc.continent,
		threshold=config["timeseries"]["threshold"],
		st_color_map=ST_COLOR_MAP
	script:
		"scripts/mlst_over_time.py"


### MRSA OVER TIME ###

rule mrsa_over_time_world:
	input:
		filtered=FILTERDIR/"all.filtered.tsv"
	output:
		table=TABLEDIR/"mrsa_over_time_world.tsv",
		plot=OUT/"{interval}"/"all"/"mrsa_over_time.all.{interval}.world.pdf"
	params:
		date_column=config["columns"]["date"],
		geography="world",
		geography_column=None,
		geography_value=None,
		continent_column=config["columns"]["continent"],
		mrsa_column=config["mrsa"]["column"],
		mrsa_value=config["mrsa"]["value"]
	script:
		"scripts/mrsa_over_time.py"


rule mrsa_over_time_region:
	input:
		filtered=FILTERDIR/"all.filtered.tsv"
	output:
		table=TABLEDIR/"mrsa_over_time_region_{region_file}.tsv",
		plot=OUT/"{interval}"/"all"/"mrsa_over_time.all.{interval}.region.{region_file}.pdf"
	params:
		date_column=config["columns"]["date"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_value=lambda wc:REGION_FILE_TO_NAME[wc.region_file],
		continent_column=config["columns"]["continent"],
		mrsa_column=config["mrsa"]["column"],
		mrsa_value=config["mrsa"]["value"]
	script:
		"scripts/mrsa_over_time.py"


rule mrsa_over_time_continent:
	input:
		filtered=FILTERDIR/"all.filtered.tsv"
	output:
		table=TABLEDIR/"mrsa_over_time_continent_{continent}.tsv",
		plot=OUT/"{interval}"/"all"/"mrsa_over_time.all.{interval}.continent.{continent}.pdf"
	params:
		date_column=config["columns"]["date"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_value=lambda wc:wc.continent,
		continent_column=config["columns"]["continent"],
		mrsa_column=config["mrsa"]["column"],
		mrsa_value=config["mrsa"]["value"]
	script:
		"scripts/mrsa_over_time.py"


### SEROTOP ###

rule serotop_world:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_serotop_world.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"serotop.{dataset}.{interval}.world.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="world",
		geography_column=None,
		geography_color_map=WORLD_COLOR_MAP,
	script:
		"scripts/serotop.py"


rule serotop_region:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_serotop_region.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"serotop.{dataset}.{interval}.region.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_color_map=REGION_COLOR_MAP,
	script:
		"scripts/serotop.py"


rule serotop_continent:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_serotop_continent.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"serotop.{dataset}.{interval}.continent.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_color_map=CONTINENT_COLOR_MAP,
	script:
		"scripts/serotop.py"


### MLST PREVALENCE PIE ###

rule mlst_prevalence_pie_world:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_prevalence_world.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_prevalence_pie.{dataset}.{interval}.world.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="world",
		geography_column=None,
		geography_value=None,
		threshold=config["pie"]["threshold"],
		st_color_map=ST_COLOR_MAP
	script:
		"scripts/mlst_prevalence_pie.py"


rule mlst_prevalence_pie_region:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_prevalence_region_{region_file}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_prevalence_pie.{dataset}.{interval}.region.{region_file}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_value=lambda wc:REGION_FILE_TO_NAME[wc.region_file],
		threshold=config["pie"]["threshold"],
		st_color_map=ST_COLOR_MAP
	script:
		"scripts/mlst_prevalence_pie.py"


rule mlst_prevalence_pie_continent:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_prevalence_continent_{continent}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_prevalence_pie.{dataset}.{interval}.continent.{continent}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_value=lambda wc:wc.continent,
		threshold=config["pie"]["threshold"],
		st_color_map=ST_COLOR_MAP
	script:
		"scripts/mlst_prevalence_pie.py"


### MLST HEATMAP ###

rule mlst_heatmap_region:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_mlst_heatmap_region.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_heatmap.{dataset}.{interval}.region.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		location_column=config["columns"]["region"],
		geography="region",
		threshold=config["heatmap"]["threshold"]
	script:
		"scripts/mlst_heatmap.py"


rule mlst_heatmap_continent:
	input:
		filtered=FILTERDIR/"{dataset}.filtered.tsv"
	output:
		table=TABLEDIR/"{dataset}_mlst_heatmap_continent.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_heatmap.{dataset}.{interval}.continent.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		location_column=config["columns"]["continent"],
		geography="continent",
		threshold=config["heatmap"]["threshold"]
	script:
		"scripts/mlst_heatmap.py"


### SEROTYPE LOGISTIC REGRESSION ###

rule serotype_logistic_regression:
	input:
		table=rules.prepare_dataset.output.filtered
	output:
		table=OUT/"{interval}"/"stats/{dataset}_serotype_logistic_regression.tsv"
	params:
		serotype_column=config["columns"]["mlst"],
		date_column=config["columns"]["date"],
		continent_column=config["columns"]["continent"],
		region_column=config["columns"]["region"],
		min_count=config["logistic_regression"]["min_count"]
	script:
		"scripts/serotype_logistic_regression.py"


rule combine_serotype_logistic_regression:
	input:
		expand(
			OUT/"{interval}"/"stats/{dataset}_serotype_logistic_regression.tsv",
			interval=YEAR_INTERVALS,
			dataset=DATASET_NAMES,
		)
	output:
		OUT/"stats/serotype_logistic_regression_all.xlsx"
	script:
		"scripts/combine_serotype_logistic_regression.py"


rule serotype_logistic_regression_figures_world:
	input:
		regression=OUT/"stats/serotype_logistic_regression_all.xlsx",
		table=rules.prepare_dataset.output.filtered
	output:
		plot=OUT/"{interval}"/"{dataset}"/"serotype_logistic_regression.{dataset}.{interval}.world.pdf"
	params:
		serotype_column=config["columns"]["mlst"],
		date_column=config["columns"]["date"],
		continent_column=config["columns"]["continent"],
		geography="world",
		geography_column=None,
		geography_value=None,
		geography_color_map=WORLD_COLOR_MAP
	script:
		"scripts/logistic_regression_figures.py"


rule serotype_logistic_regression_figures_region:
	input:
		regression=OUT/"stats/serotype_logistic_regression_all.xlsx",
		table=rules.prepare_dataset.output.filtered
	output:
		plot=OUT/"{interval}"/"{dataset}"/"serotype_logistic_regression.{dataset}.{interval}.region.{region_file}.pdf"
	params:
		serotype_column=config["columns"]["mlst"],
		date_column=config["columns"]["date"],
		continent_column=config["columns"]["continent"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_value=lambda wc:REGION_FILE_TO_NAME[wc.region_file],
		geography_color_map=REGION_COLOR_MAP
	script:
		"scripts/logistic_regression_figures.py"


rule serotype_logistic_regression_figures_continent:
	input:
		regression=OUT/"stats/serotype_logistic_regression_all.xlsx",
		table=rules.prepare_dataset.output.filtered
	output:
		plot=OUT/"{interval}"/"{dataset}"/"serotype_logistic_regression.{dataset}.{interval}.continent.{continent}.pdf"
	params:
		serotype_column=config["columns"]["mlst"],
		date_column=config["columns"]["date"],
		continent_column=config["columns"]["continent"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_value=lambda wc:wc.continent,
		geography_color_map=CONTINENT_COLOR_MAP
	script:
		"scripts/logistic_regression_figures.py"


### TAR OVER TIME ###

rule tar_over_time_world:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_tar_over_time_world.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"tar_over_time.{dataset}.{interval}.world.pdf"
	params:
		date_column=config["columns"]["date"],
		geography="world",
		geography_column=None,
		geography_value=None,
		tar_colors=TAR_COLORS
	script:
		"scripts/tar_over_time.py"


rule tar_over_time_region:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_tar_over_time_region_{region_file}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"tar_over_time.{dataset}.{interval}.region.{region_file}.pdf"
	params:
		date_column=config["columns"]["date"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_value=lambda wc:REGION_FILE_TO_NAME[wc.region_file],
		tar_colors=TAR_COLORS
	script:
		"scripts/tar_over_time.py"


rule tar_over_time_continent:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_tar_over_time_continent_{continent}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"tar_over_time.{dataset}.{interval}.continent.{continent}.pdf"
	params:
		date_column=config["columns"]["date"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_value=lambda wc:wc.continent,
		tar_colors=TAR_COLORS
	script:
		"scripts/tar_over_time.py"


### MLST TAR PREVALENCE PIE ###

rule mlst_tar_prevalence_pie_world:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_tar_prevalence_world.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_tar_prevalence_pie.{dataset}.{interval}.world.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="world",
		geography_column=None,
		geography_value=None,
		threshold=config["pie"]["threshold"],
		st_color_map=ST_COLOR_MAP,
		tar_colors=TAR_COLORS
	script:
		"scripts/mlst_tar_prevalence_pie.py"


rule mlst_tar_prevalence_pie_region:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_tar_prevalence_region_{region_file}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_tar_prevalence_pie.{dataset}.{interval}.region.{region_file}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_value=lambda wc:REGION_FILE_TO_NAME[wc.region_file],
		threshold=config["pie"]["threshold"],
		st_color_map=ST_COLOR_MAP,
		tar_colors=TAR_COLORS
	script:
		"scripts/mlst_tar_prevalence_pie.py"


rule mlst_tar_prevalence_pie_continent:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_tar_prevalence_continent_{continent}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_tar_prevalence_pie.{dataset}.{interval}.continent.{continent}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_value=lambda wc:wc.continent,
		threshold=config["pie"]["threshold"],
		st_color_map=ST_COLOR_MAP,
		tar_colors=TAR_COLORS
	script:
		"scripts/mlst_tar_prevalence_pie.py"


### MLST TAR OVER TIME ###

rule mlst_tar_over_time_world:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_tar_over_time_world.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_tar_over_time.{dataset}.{interval}.world.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		date_column=config["columns"]["date"],
		geography="world",
		geography_column=None,
		geography_value=None,
		threshold=config["pie"]["threshold"],
		mlst_tar_colors=MLST_TAR_COLORS
	script:
		"scripts/mlst_tar_over_time.py"


rule mlst_tar_over_time_region:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_tar_over_time_region_{region_file}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_tar_over_time.{dataset}.{interval}.region.{region_file}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		date_column=config["columns"]["date"],
		geography="region",
		geography_column=config["columns"]["region"],
		geography_value=lambda wc:REGION_FILE_TO_NAME[wc.region_file],
		threshold=config["pie"]["threshold"],
		mlst_tar_colors=MLST_TAR_COLORS
	script:
		"scripts/mlst_tar_over_time.py"


rule mlst_tar_over_time_continent:
	input:
		filtered=rules.prepare_dataset.output.filtered
	output:
		table=TABLEDIR/"{dataset}_mlst_tar_over_time_continent_{continent}.tsv",
		plot=OUT/"{interval}"/"{dataset}"/"mlst_tar_over_time.{dataset}.{interval}.continent.{continent}.pdf"
	params:
		mlst_column=config["columns"]["mlst"],
		date_column=config["columns"]["date"],
		geography="continent",
		geography_column=config["columns"]["continent"],
		geography_value=lambda wc:wc.continent,
		threshold=config["pie"]["threshold"],
		mlst_tar_colors=MLST_TAR_COLORS
	script:
		"scripts/mlst_tar_over_time.py"