import pandas as pd

from utils import load_data, get_mrsa


def summarize(df, geography, name, mlst_column, mrsa_column, mrsa_value):

	mrsa = get_mrsa(
		df,
		mrsa_column,
		mrsa_value,
	)

	return {
		"geography": geography,
		"name": name,
		"genomes": len(df),
		"MLSTs": df[mlst_column].dropna().nunique(),
		"MRSA": int(mrsa.sum()),
		"non_MRSA": int((~mrsa).sum()),
	}


def run_analysis(smk):

	df = load_data(
		smk.input.filtered,
		smk.params.mlst_column
	)

	mlst_column = smk.params.mlst_column

	results = []

	# world
	results.append(
		summarize(
			df,
			"world",
			"world",
			mlst_column,
			smk.params.mrsa_column,
			smk.params.mrsa_value,
		)
	)

	# continents
	for continent, sub in df.groupby(
		smk.params.continent_column
	):

		results.append(
			summarize(
				sub,
				"continent",
				str(continent),
				mlst_column,
				smk.params.mrsa_column,
				smk.params.mrsa_value,
			)
		)

	# regions
	for region, sub in df.groupby(
		smk.params.region_column
	):

		results.append(
			summarize(
				sub,
				"region",
				str(region),
				mlst_column,
				smk.params.mrsa_column,
				smk.params.mrsa_value,
			)
		)

	pd.DataFrame(results).to_csv(
		smk.output.table,
		sep="\t",
		index=False,
	)


if __name__ == "__main__":
	run_analysis(snakemake)  # type: ignore