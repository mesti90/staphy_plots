from pathlib import Path
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import pandas as pd
from utils import load_data


CATEGORIES=["none","tarM","tarM+tarP","tarM+tarP+tarS","tarM+tarS","tarP","tarP+tarS","tarS"]


def write_empty(outfile,table,title,columns):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(10,6))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)
	pd.DataFrame(columns=columns).to_csv(table,sep="\t",index=False)


def prepare_tar(df):
	tar_columns={category:f"tar|{category}" for category in CATEGORIES}
	missing=[column for column in tar_columns.values() if column not in df.columns]
	if missing:
		raise ValueError("Missing tar columns: "+", ".join(missing))
	tar=df[list(tar_columns.values())].apply(pd.to_numeric,errors="coerce").fillna(0)
	tar.columns=CATEGORIES
	has_tar=tar.sum(axis=1)>0
	df=df.loc[has_tar].copy()
	tar=tar.loc[has_tar]
	df["tar_combination"]=tar.idxmax(axis=1)
	return df


def plot_tar(df,outfile,title,tar_colors):
	counts=df.groupby(["year","tar_combination"]).size().reset_index(name="count")
	totals=counts.groupby("year")["count"].sum()
	counts["prevalence"]=counts["count"]/counts["year"].map(totals)*100
	wide=(
		counts.pivot(index="year",columns="tar_combination",values="prevalence")
		.reindex(columns=CATEGORIES)
		.fillna(0)
	)
	plot_data=wide.loc[:,(wide>0).any(axis=0)]
	fig,ax=plt.subplots(figsize=(10,6))
	colors=[tar_colors[category] for category in plot_data.columns]
	ax.stackplot(
		plot_data.index,
		*[plot_data[category] for category in plot_data.columns],
		labels=plot_data.columns,
		colors=colors,
	)
	ax.set_xlabel("Year")
	ax.set_ylabel("Prevalence (%)")
	ax.set_ylim(0,100)
	ax.set_xticks(sorted(plot_data.index)[::2])
	ax.set_title(title)
	ax.legend(title="Tar combination",loc="center left",bbox_to_anchor=(1,0.5))
	fig.tight_layout()
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close(fig)
	return wide


def run_analysis(smk):
	df=load_data(smk.input.filtered)
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
			["year"]+CATEGORIES,
		)
		return
	df=prepare_tar(df)
	if df.empty:
		if geography=="world":
			raise ValueError("No tar-positive data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No tar-positive data for {geography}: {geography_value}",
			["year"]+CATEGORIES,
		)
		return
	df["year"]=pd.to_datetime(df[smk.params.date_column],errors="coerce").dt.year
	df=df.dropna(subset=["year","tar_combination"])
	df["year"]=df["year"].astype(int)
	if df.empty:
		if geography=="world":
			raise ValueError("No valid temporal data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No valid temporal data for {geography}: {geography_value}",
			["year"]+CATEGORIES,
		)
		return
	wide=plot_tar(
		df,
		smk.output.plot,
		f"Tar combinations over time - {geography}"
		+("" if geography=="world" else f": {geography_value}")
		+f" ({smk.wildcards.dataset}, {smk.wildcards.interval})",
		smk.params.tar_colors,
	)
	wide.to_csv(smk.output.table,sep="\t")


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore