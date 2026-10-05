from pathlib import Path

import pandas as pd
from openpyxl.utils import get_column_letter


def run_analysis(smk):

	results = []

	for file in smk.input:

		path = Path(file)

		# parse:
		# .../{interval}/stats/{dataset}_serotype_logistic_regression.tsv

		interval = path.parents[1].name

		dataset = ( path.stem .replace("_serotype_logistic_regression", "") )

		df = pd.read_csv( path, sep="\t", )
		df["dataset"] = dataset
		df["interval"] = interval

		results.append(df)

	out = pd.concat( results, ignore_index=True, )
	# optional: nicer ordering
	first = [ "dataset", "interval", "serotype", "scope", "count", "n_total", "predictor", "estimate", "odds_ratio", "ci_lower", "ci_upper", "p_value", "p_value_corr", ]

	out["trend_interpretation"] = out.apply( interpret_trend, axis=1, )
	out = out[ [c for c in first if c in out.columns] + [ c for c in out.columns if c not in first ] ]
	out = out.sort_values(["serotype", "dataset", "interval", "predictor"])
	Path(smk.output[0]).parent.mkdir( parents=True, exist_ok=True, )
	legend = pd.DataFrame(
		[
			{
				"Trend interpretation": "significant increasing trend",
				"Meaning": "Odds ratio ≥ 1.05 and FDR-adjusted p-value < 0.05. The serotype became significantly more prevalent over time."
			},
			{
				"Trend interpretation": "significant decreasing trend",
				"Meaning": "Odds ratio ≤ 0.95 and FDR-adjusted p-value < 0.05. The serotype became significantly less prevalent over time."
			},
			{
				"Trend interpretation": "significant weak temporal trend",
				"Meaning": "FDR-adjusted p-value < 0.05, but the annual effect size is small (0.95 < OR < 1.05)."
			},
			{
				"Trend interpretation": "non-significant increasing trend",
				"Meaning": "Odds ratio > 1.05, but the trend is not statistically significant after FDR correction."
			},
			{
				"Trend interpretation": "non-significant decreasing trend",
				"Meaning": "Odds ratio < 0.95, but the trend is not statistically significant after FDR correction."
			},
			{
				"Trend interpretation": "non-significant weak temporal trend",
				"Meaning": "No statistically significant temporal trend and only a small estimated effect (0.95 ≤ OR ≤ 1.05)."
			},
		]
	)

	with pd.ExcelWriter(smk.output[0], engine="openpyxl") as writer:
		out.to_excel(writer, sheet_name="Results", index=False)
		legend.to_excel(writer, sheet_name="Trend interpretation", index=False)

		for worksheet in writer.sheets.values():
			for column_cells in worksheet.columns:
				max_length = max(
					len(str(cell.value)) if cell.value is not None else 0
					for cell in column_cells
				)
				column_letter = get_column_letter(column_cells[0].column)
				worksheet.column_dimensions[column_letter].width = min(max_length + 2, 80)

def interpret_trend(row):

	or_value = row["odds_ratio"]
	p = row["p_value_corr"]

	if p >= 0.05:
		if or_value > 1.05:
			return "non-significant increasing trend"
		elif or_value < 0.95:
			return "non-significant decreasing trend"
		else:
			return "non-significant weak temporal trend"

	if or_value >= 1.05:
		return "significant increasing trend"

	if or_value <= 0.95:
		return "significant decreasing trend"

	return "significant weak temporal trend"

if __name__ == "__main__":
	run_analysis(snakemake) # type: ignore