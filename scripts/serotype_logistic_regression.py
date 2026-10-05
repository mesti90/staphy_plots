from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf # type: ignore
from statsmodels.stats.multitest import multipletests # type: ignore

from utils import load_data


def fit_model(df, formula, predictor, model_name, serotype_column, min_count):
	results = []
	for serotype, sub in df.groupby(serotype_column):
		count = len(sub)
		if count < min_count:
			continue
		work = df.copy()
		work["present"] = ( work[serotype_column] == serotype ).astype(int)
		try:
			model = smf.logit( formula=formula, data=work, ).fit(disp=False)
			null_model = smf.logit( formula="present ~ 1", data=work, ).fit(disp=False)
			mcfadden_r2 = 1 - (model.llf / null_model.llf)
		except Exception:
			continue
		if predictor not in model.params.index:
			continue
		ci = model.conf_int().loc[predictor]
		results.append({
			"serotype": serotype,
			"scope": "world",
			"count": count,
			"n_total": len(work),
			"predictor": model_name,
			#"r_squared": mcfadden_r2,
			"estimate": model.params[predictor],
			"odds_ratio": np.exp(model.params[predictor]),
			"ci_lower": np.exp(ci[0]),
			"ci_upper": np.exp(ci[1]),
			"p_value": model.pvalues[predictor],
		})

	results = pd.DataFrame(results)
	if len(results):
		results["p_value_corr"] = multipletests( results["p_value"], method="fdr_bh", )[1]
	return results


def run_analysis(smk):

	df = load_data( smk.input.table, smk.params.serotype_column, )
	df = df.dropna( subset=[ smk.params.serotype_column, smk.params.date_column, smk.params.continent_column, smk.params.region_column, ] ).copy()
	df["collection_year"] = ( df[smk.params.date_column] .astype(str) .str.extract(r"^(\d{4})")[0] .astype(float) )
	models = [
		( "present ~ collection_year", "collection_year", "collection_year", ),
		( "present ~ collection_year + C(continent)", "collection_year", "collection_year, continent", ),
		( "present ~ collection_year + C(UNRegion)", "collection_year", "collection_year, UNRegion", ),
	]
	all_results = []
	for formula, predictor, model_name in models:
		work = df.rename( columns={ smk.params.continent_column: "continent", smk.params.region_column: "UNRegion", } )
		res = fit_model( work, formula, predictor, model_name, smk.params.serotype_column, smk.params.min_count, )
		all_results.append(res)
	results = pd.concat( all_results, ignore_index=True, )
	Path(smk.output[0]).parent.mkdir( parents=True, exist_ok=True, )
	results.to_csv( smk.output[0], sep="\t", index=False, )




if __name__ == "__main__":
	run_analysis(snakemake)  # type: ignore