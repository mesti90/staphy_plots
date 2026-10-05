from pathlib import Path
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import pandas as pd
from utils import load_data


def write_empty(outfile,table,title,columns):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(8,8))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)
	pd.DataFrame(columns=columns).to_csv(table,sep="\t",index=False)


def prepare_pie(df,mlst_column,threshold):
	df=df.dropna(subset=[mlst_column]).copy()
	counts=df[mlst_column].value_counts()
	all_types=len(counts)
	cum=counts.cumsum()/counts.sum()
	keep=list(cum[cum<=threshold].index)
	if not keep:
		keep=[counts.index[0]]
	out=df[mlst_column].where(df[mlst_column].isin(keep),"Other")
	out=out.value_counts().reset_index()
	out.columns=["mlst","count"]
	if "Other" in out["mlst"].values:
		other_types=all_types-len(keep)
		out.loc[out["mlst"]=="Other","mlst"]=f"Other ({other_types})"
	out["percentage"]=out["count"]/out["count"].sum()*100
	out["is_other"]=out["mlst"].str.startswith("Other")
	return out.sort_values(["is_other","count"],ascending=[True,False]).drop(columns="is_other")


def plot_pie(df,outfile,title,st_color_map):
	df=df.copy()
	df["is_other"]=df["mlst"].str.startswith("Other")
	df=df.sort_values(["is_other","count"],ascending=[True,False]).drop(columns="is_other")
	colors=[
		"#808080" if str(mlst).startswith("Other") else st_color_map.get(str(mlst),"#808080")
		for mlst in df["mlst"]
	]
	plt.figure(figsize=(8,8))
	plt.pie(df["count"],labels=None,colors=colors,startangle=90,wedgeprops={"edgecolor":"white"})
	plt.legend(
		labels=[f"{mlst} (n={count}, {perc:.1f}%)" for mlst,count,perc in zip(df["mlst"],df["count"],df["percentage"])],
		title=f"MLST (Total: n={df['count'].sum()})",
		bbox_to_anchor=(1.05,1),loc="upper left",fontsize=8,
	)
	plt.title(title)
	plt.tight_layout()
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	plt.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close()


def run_analysis(smk):
	df=load_data(smk.input.filtered,smk.params.mlst_column)
	geography=smk.params.geography
	geography_value=smk.params.geography_value
	if smk.params.geography_column is not None:
		df=df[df[smk.params.geography_column]==geography_value].copy()
	if df.empty:
		if geography=="world":
			raise ValueError("No data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No data for {geography}: {geography_value}",
			["mlst","count","percentage"],
		)
		return
	pie=prepare_pie(df,smk.params.mlst_column,float(smk.params.threshold))
	if pie.empty:
		if geography=="world":
			raise ValueError("No MLST data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No MLST data for {geography}: {geography_value}",
			["mlst","count","percentage"],
		)
		return
	pie.to_csv(smk.output.table,sep="\t",index=False)
	title=f"MLST prevalence - {geography}"
	if geography!="world":
		title+=f": {geography_value}"
	title+=f" - genomes: {smk.wildcards.dataset} ({smk.wildcards.interval})"
	plot_pie(pie,smk.output.plot,title,smk.params.st_color_map)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore