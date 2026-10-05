from pathlib import Path
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import pandas as pd
from utils import load_data,get_mrsa


def write_empty(outfile,table,title,columns):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(10,6))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)
	pd.DataFrame(columns=columns).to_csv(table,sep="\t",index=False)


def mrsa_frequency(df,date_column,mrsa_column,mrsa_value,geography_column=None,geography_type="world",geography_value=None):
	df=df.copy()
	df["year"]=pd.to_numeric(df[date_column].astype(str).str[:4],errors="coerce")
	df["MRSA"]=get_mrsa(df,mrsa_column,mrsa_value)
	if geography_column is None:
		df["geography"]="world"
	else:
		df=df[df[geography_column]==geography_value].copy()
		df["geography"]=geography_value
	df["geography_type"]=geography_type
	return df.groupby(["geography_type","geography","year","MRSA"]).size().reset_index(name="count")


def plot_mrsa(df,outfile,title):
	p=df.groupby(["year","MRSA"])["count"].sum().unstack(fill_value=0).sort_index()
	years=p.index.to_numpy()
	x=range(len(years))
	width=0.4
	plt.figure(figsize=(10,6))
	if False in p.columns:
		plt.bar([i-width/2 for i in x],p[False],width=width,label="non-MRSA")
	if True in p.columns:
		plt.bar([i+width/2 for i in x],p[True],width=width,label="MRSA")
	plt.xticks(list(x),years.astype(int),rotation=45)
	plt.xlabel("Year")
	plt.ylabel("# genomes")
	plt.title(title,pad=15)
	plt.grid(True,axis="y",linestyle="--",linewidth=0.5,alpha=0.4)
	plt.legend()
	plt.tight_layout()
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	plt.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close()


def run_analysis(smk):
	df=load_data(smk.input.filtered)
	geography=smk.params.geography
	geography_value=smk.params.geography_value
	freq=mrsa_frequency(
		df,
		smk.params.date_column,
		smk.params.mrsa_column,
		smk.params.mrsa_value,
		smk.params.geography_column,
		geography,
		geography_value,
	)
	if freq.empty:
		if geography=="world":
			raise ValueError("No data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No data for {geography}: {geography_value}",
			["geography_type","geography","year","MRSA","count"],
		)
		return
	freq.to_csv(smk.output.table,sep="\t",index=False)
	title=f"MRSA prevalence over time - {geography}"
	if geography!="world":
		title+=f": {geography_value}"
	title+=f" - {smk.wildcards.interval}"
	plot_mrsa(freq,smk.output.plot,title)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore