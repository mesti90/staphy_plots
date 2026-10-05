from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
	snakemake=None

import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore
import pandas as pd
from matplotlib.colors import PowerNorm,LinearSegmentedColormap # type: ignore
from utils import load_data


def prepare_matrix(df,location_column,mlst_column,threshold):
	df=df.dropna(subset=[mlst_column,location_column]).copy()
	counts=(
		df.groupby([location_column,mlst_column])
		.size()
		.reset_index(name="count")
	)
	total=counts.groupby(location_column)["count"].transform("sum")
	counts["percent"]=counts["count"]/total*100
	global_counts=(
		counts.groupby(mlst_column)["count"]
		.sum()
		.sort_values(ascending=False)
	)
	cum=global_counts.cumsum()/global_counts.sum()
	keep=list(cum[cum<=threshold].index)
	if not keep:
		keep=[global_counts.index[0]]
	counts["mlst_plot"]=counts[mlst_column].where(
		counts[mlst_column].isin(keep),
		"Other"
	)
	counts=(
		counts.groupby([location_column,"mlst_plot"])["percent"]
		.sum()
		.reset_index()
	)
	mat=(
		counts
		.pivot(
			index=location_column,
			columns="mlst_plot",
			values="percent"
		)
		.fillna(0)
	)
	if "Other" in mat.columns:
		mat=mat[
			[x for x in mat.columns if x!="Other"]+["Other"]
		]
	return mat


def plot_heatmap(mat,outfile,title):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	if len(mat)<2:
		fig,ax=plt.subplots(figsize=(10,6))
		ax.axis("off")
		ax.text(
			0.5,
			0.5,
			f"Not enough locations for heatmap\n{title}",
			ha="center",
			va="center",
		)
		fig.savefig(outfile,bbox_inches="tight")
		plt.close(fig)
		return
	cmap=LinearSegmentedColormap.from_list(
		"grey_blue",
		[
			(0.0,"lightgrey"),
			(0.001,"lightgrey"),
			(0.01,"white"),
			(0.1,"lightblue"),
			(0.5,"royalblue"),
			(1.0,"darkblue"),
		]
	)
	g=sns.clustermap(
		mat,
		method="average",
		metric="braycurtis",
		row_cluster=True,
		col_cluster=False,
		cmap=cmap,
		norm=PowerNorm(gamma=0.4),
		vmin=0,
		vmax=100,
		linewidths=0.2,
		figsize=(14,max(5,len(mat)*0.35)),
		xticklabels=True,
		yticklabels=True,
		cbar_kws={"label":"% genomes"},
	)
	g.ax_heatmap.tick_params(
		axis="y",
		labelsize=8,
		pad=8,
	)
	g.ax_heatmap.tick_params(
		axis="x",
		labelrotation=90,
		labelsize=8,
	)
	g.ax_heatmap.set_title(title,pad=20)
	g.ax_heatmap.set_xlabel("MLST")
	g.ax_heatmap.set_ylabel("")
	g.savefig(
		outfile,
		dpi=300,
		bbox_inches="tight",
	)
	plt.close(g.fig)


def run_analysis(smk):
	df=load_data(
		smk.input.filtered,
		smk.params.mlst_column
	)
	location_column=smk.params.location_column
	mat=prepare_matrix(
		df,
		location_column,
		smk.params.mlst_column,
		float(smk.params.threshold),
	)
	mat.to_csv(
		smk.output.table,
		sep="\t",
	)
	plot_heatmap(
		mat,
		smk.output.plot,
		f"MLST distribution by {smk.params.geography} - "
		f"genomes: {smk.wildcards.dataset} "
		f"({smk.wildcards.interval})",
	)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore