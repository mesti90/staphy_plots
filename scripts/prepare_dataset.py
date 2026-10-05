from pathlib import Path

from utils import load_data, filter_year


def run_analysis(smk):

	df = load_data(
		smk.input.tsv,
		"mlst"
	)

	df = filter_year(
		df,
		int(smk.params.minyear),
		int(smk.params.maxyear),
		smk.params.date_column,
	)

	# Remove genomes without valid dates or UNRegion
	df = df.dropna(
		subset=[
			"year",
			smk.params.region_column,
		]
	)

	# Add continent information from UNRegion
	continent_column = smk.params.continent_column

	df[continent_column] = df[
		smk.params.region_column
	].map(
		smk.params.region_to_continent
	)

	# Check for unmapped regions
	missing = df.loc[
		df[continent_column].isna(),
		smk.params.region_column
	].unique()

	if len(missing):
		raise ValueError(
			f"Missing continent mapping for regions: {list(missing)}"
		)

	Path(smk.output.filtered).parent.mkdir(
		parents=True,
		exist_ok=True,
	)

	df.to_csv(
		smk.output.filtered,
		sep="\t",
		index=False,
	)


if __name__=="__main__":
	run_analysis(snakemake)  # type: ignore