from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
	snakemake=None

import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import pandas as pd
from utils import load_data


def calculate_saturation(df,geography_column,mlst_column):
	results=[]
	if geography_column is None:
		groups=[("world",df)]
	else:
		groups=df.groupby(geography_column)
	for geography,sub in groups:
		sub=sub.dropna(subset=[mlst_column])
		if len(sub)==0:
			continue
		results.append({
			"geography":geography,
			"mlst":None,
			"count":0,
			"types":0,
			"quantile":0.0
		})
		counts=sub[mlst_column].value_counts()
		total=counts.sum()
		cumulative=0
		for order,(mlst,count) in enumerate(counts.items(),start=1):
			cumulative+=count
			results.append({
				"geography":geography,
				"mlst":mlst,
				"count":count,
				"types":order,
				"quantile":cumulative/total
			})
	return pd.DataFrame(results)


def plot_saturation(df,outfile,title,geography_color_map):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(10,6))
	if df.empty:
		ax.axis("off")
		ax.text(0.5,0.5,"No data",ha="center",va="center")
		fig.savefig(outfile,bbox_inches="tight")
		plt.close(fig)
		return
	for geography in df["geography"].unique():
		sub=df[df["geography"]==geography].sort_values("quantile")
		color=geography_color_map.get(str(geography),"#808080")
		ax.step(
			sub["quantile"],
			sub["types"],
			where="post",
			label=str(geography),
			color=color,
		)
	ax.set_xlabel("Fraction of genomes")
	ax.set_ylabel("# MLST types")
	ax.set_title(title)
	ax.set_xlim(0,1)
	ax.set_ylim(0,None)
	ax.grid(True,linestyle="--",linewidth=0.5,alpha=0.4)
	ax.legend(
		fontsize=8,
		bbox_to_anchor=(1.05,1),
		loc="upper left"
	)
	plt.tight_layout()
	fig.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close(fig)


def run_analysis(smk):
	df=load_data(smk.input.filtered)
	geography_column=smk.params.geography_column
	if geography_column is None:
		df["world"]="world"
		geography_column="world"
	saturation=calculate_saturation(
		df,
		geography_column,
		smk.params.mlst_column,
	)
	saturation.to_csv(
		smk.output.table,
		sep="\t",
		index=False
	)
	plot_saturation(
		saturation,
		smk.output.plot,
		f"MLST accumulation - {smk.params.geography} - "
		f"genomes: {smk.wildcards.dataset} "
		f"({smk.wildcards.interval})",
		smk.params.geography_color_map,
	)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore