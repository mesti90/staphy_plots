from pathlib import Path
import math
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
from matplotlib.backends.backend_pdf import PdfPages # type: ignore
import pandas as pd
from utils import load_data


def write_empty(outfile,title):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(8,2))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)


def interpret_arrow(estimate):
	if estimate>0:
		return "↑","red"
	if estimate<0:
		return "↓","forestgreen"
	return "","black"


def yearly_prevalence(df,serotype_column,serotype):
	work=df.assign(present=df[serotype_column]==serotype)
	stats=(
		work.groupby("collection_year")
		.agg(n_total=("present","size"),n_present=("present","sum"))
		.reset_index()
	)
	stats["prevalence"]=stats["n_present"]/stats["n_total"]
	return stats


def prepare_plot_data(df,serotype,serotype_column):
	if len(df)==0:
		return None
	first=df.loc[df[serotype_column]==serotype,"collection_year"]
	if len(first)==0:
		return None
	first_year=first.min()
	df=df[df["collection_year"]>=first_year].copy()
	stats=yearly_prevalence(df,serotype_column,serotype)
	stats=stats[stats["n_present"]>0].copy()
	return stats if len(stats) else None


def draw_panel(ax,stats,serotype,p_value,estimate,geography,geography_value):
	arrow,color=interpret_arrow(estimate)
	ax.bar(stats["collection_year"].astype(str),stats["prevalence"],color="darkgrey")
	ymax=max(stats["prevalence"].max(),1)
	for _,row in stats.iterrows():
		ax.text(str(int(row["collection_year"])),row["prevalence"]+ymax*0.02,str(int(row["n_present"])),rotation=90,ha="center",va="bottom",fontsize=7)
	title=f"{serotype} ({geography}"
	if geography!="world":
		title+=f": {geography_value}"
	title+=f")\nFDR={p_value:.3g}"
	ax.set_title(title,fontsize=9)
	if arrow:
		ax.text(0.98,1.02,arrow,transform=ax.transAxes,ha="right",va="bottom",color=color,fontsize=14,fontweight="bold")
	ax.set_xlabel("")
	ax.set_ylabel("Relative prevalence")
	ax.tick_params(axis="x",labelrotation=90,labelsize=8)
	ax.tick_params(axis="y",labelsize=8)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.set_ylim(0,ymax*1.15)
	ax.grid(axis="y",alpha=0.3)


def run_analysis(smk):
	df=load_data(smk.input.table,smk.params.serotype_column)
	geography=smk.params.geography
	geography_value=smk.params.geography_value
	if smk.params.geography_column is not None:
		df=df[df[smk.params.geography_column]==geography_value].copy()
	if df.empty:
		if geography=="world":
			raise ValueError("No data for world")
		write_empty(smk.output.plot,f"No data for {geography}: {geography_value}")
		return
	df["collection_year"]=df[smk.params.date_column].astype(str).str.extract(r"^(\d{4})")[0]
	df=df.dropna(subset=["collection_year"]).copy()
	df["collection_year"]=df["collection_year"].astype(int)
	if df.empty:
		if geography=="world":
			raise ValueError("No valid collection years for world")
		write_empty(smk.output.plot,f"No valid collection years for {geography}: {geography_value}")
		return
	results=pd.read_excel(smk.input.regression,sheet_name="Results")
	results=results[
		(results["dataset"]==smk.wildcards.dataset)
		&(results["interval"]==smk.wildcards.interval)
		&(results["predictor"]=="collection_year")
		&(results["scope"]==geography)
		&(results["p_value_corr"]<0.05)
	].copy()
	Path(smk.output.plot).parent.mkdir(parents=True,exist_ok=True)
	if len(results)==0:
		write_empty(smk.output.plot,"No significant temporal trends.")
		return
	results["_st_num"]=pd.to_numeric(results["serotype"].str.extract(r"(\d+)")[0],errors="coerce")
	results=results.sort_values(["_st_num","serotype"],na_position="last").drop(columns="_st_num")
	with PdfPages(smk.output.plot) as pdf:
		for start in range(0,len(results),12):
			page=results.iloc[start:start+12]
			if len(page)==4:
				ncols=2
			elif len(page)<=3:
				ncols=max(1,len(page))
			else:
				ncols=3
			nrows=math.ceil(len(page)/ncols)
			fig,axes=plt.subplots(nrows,ncols,figsize=(6*ncols,3.5*nrows),squeeze=False)
			axes=axes.ravel()
			for ax,(_,row) in zip(axes,page.iterrows()):
				stats=prepare_plot_data(df,row["serotype"],smk.params.serotype_column)
				if stats is None:
					ax.axis("off")
					continue
				draw_panel(ax,stats,row["serotype"],row["p_value_corr"],row["estimate"],geography,geography_value)
			for ax in axes[len(page):]:
				ax.axis("off")
			title=f"Significant temporal MLST trends\n{smk.wildcards.dataset} | {smk.wildcards.interval}"
			if geography!="world":
				title+=f" | {geography}: {geography_value}"
			fig.suptitle(title,fontsize=14)
			fig.supxlabel("Collection year")
			fig.supylabel("Relative prevalence")
			fig.tight_layout(rect=(0.03,0.03,1,0.95))
			pdf.savefig(fig)
			plt.close(fig)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore