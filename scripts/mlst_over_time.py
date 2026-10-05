from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
	snakemake=None

import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import matplotlib.ticker as ticker # type: ignore
import pandas as pd
from utils import load_data


def write_empty(outfile,table,title,columns):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(10,6))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)
	pd.DataFrame(columns=columns).to_csv(table,sep="\t",index=False)


def mlst_frequency(df,mlst_column,geography_column=None,geography_type="world",geography_value=None):
	df=df.copy()
	if geography_column is None:
		df["geography"]="world"
	else:
		df=df[df[geography_column]==geography_value].copy()
		df["geography"]=geography_value
	df["geography_type"]=geography_type
	out=(
		df.groupby(["geography_type","geography","year",mlst_column])
		.size()
		.reset_index(name="count")
	)
	out["ratio"]=(
		out.groupby(["geography_type","geography","year"])["count"]
		.transform(lambda x:x/x.sum())
	)
	return out


def select_common_types(freq,mlst_column,threshold):
	global_freq=freq.groupby(mlst_column)["count"].sum().sort_values(ascending=False)
	cum=global_freq.cumsum()/global_freq.sum()
	keep=list(cum[cum<=threshold].index)
	if not keep:
		keep=[global_freq.index[0]]
	return keep


def plot_area(freq,keep,mlst_column,st_color_map,outfile,title):
	df=freq.copy()
	df["mlst_plot"]=df[mlst_column].where(df[mlst_column].isin(keep),"Other")
	plot=df.groupby(["year","mlst_plot"])["ratio"].sum().unstack(fill_value=0).sort_index()

	def mlst_order(x):
		if x=="Other":
			return float("inf")
		try:
			return int(str(x).replace("ST",""))
		except ValueError:
			return float("inf")

	plot=plot[sorted(plot.columns,key=mlst_order)]
	colors=[
		"#808080" if st=="Other" else st_color_map.get(str(st),"#808080")
		for st in plot.columns
	]
	plt.figure(figsize=(10,6))
	plt.stackplot(plot.index,plot.T,labels=plot.columns,colors=colors)
	plt.xlabel("Year")
	plt.ylabel("Relative prevalence")
	plt.title(title,pad=15)
	ax=plt.gca()
	ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
	handles,labels=ax.get_legend_handles_labels()
	ax.legend(handles,labels,bbox_to_anchor=(1.02,1),loc="upper left",fontsize=8)
	plt.tight_layout()
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	plt.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close()


def run_analysis(smk):
	df=load_data(smk.input.filtered,smk.params.mlst_column)
	mlst_column=smk.params.mlst_column
	geography=smk.params.geography
	geography_column=smk.params.geography_column
	geography_value=smk.params.geography_value
	freq=mlst_frequency(df,mlst_column,geography_column,geography,geography_value)
	if freq.empty:
		if geography=="world":
			raise ValueError("No data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No data for {geography}: {geography_value}",
			["geography_type","geography","year",mlst_column,"count","ratio"],
		)
		return
	keep=select_common_types(freq,mlst_column,float(smk.params.threshold))
	freq.to_csv(smk.output.table,sep="\t",index=False)
	plot_area(
		freq,
		keep,
		mlst_column,
		smk.params.st_color_map,
		smk.output.plot,
		f"MLST prevalence over time - {geography}"
		+("" if geography=="world" else f": {geography_value}")
		+f" - {smk.wildcards.dataset} ({smk.wildcards.interval})",
	)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore