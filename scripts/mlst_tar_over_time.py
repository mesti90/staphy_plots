from pathlib import Path
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import pandas as pd
from utils import load_data


TAR_CATEGORIES=["none","tarM","tarM+tarP","tarM+tarP+tarS","tarM+tarS","tarP","tarP+tarS","tarS"]


def write_empty(outfile,table,title,columns):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(11,7))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)
	pd.DataFrame(columns=columns).to_csv(table,sep="\t",index=False)


def prepare_tar(df):
	tar_columns={category:f"tar|{category}" for category in TAR_CATEGORIES}
	missing=[column for column in tar_columns.values() if column not in df.columns]
	if missing:
		raise ValueError("Missing tar columns: "+", ".join(missing))
	tar=df[list(tar_columns.values())].apply(pd.to_numeric,errors="coerce").fillna(0)
	tar.columns=TAR_CATEGORIES
	if (tar.sum(axis=1)>1).any():
		raise ValueError("Some rows have multiple tar combination columns set to 1")
	has_tar=tar.sum(axis=1)>0
	df=df.loc[has_tar].copy()
	tar=tar.loc[has_tar]
	df["tar_combination"]=tar.idxmax(axis=1)
	return df


def prepare_data(df,mlst_column,date_column,threshold):
	df=df.dropna(subset=[mlst_column,date_column]).copy()
	df["year"]=pd.to_datetime(df[date_column],errors="coerce").dt.year
	df=df.dropna(subset=["year"])
	df["year"]=df["year"].astype(int)
	df["mlst_tar"]=df[mlst_column].astype(str)+" + "+df["tar_combination"].astype(str)
	counts=df["mlst_tar"].value_counts()
	cumulative=counts.cumsum()/counts.sum()
	keep=list(cumulative[cumulative<=threshold].index)
	if not keep:
		keep=[counts.index[0]]
	df["mlst_tar"]=df["mlst_tar"].where(df["mlst_tar"].isin(keep),"Other")
	return df,keep


def plot_over_time(df,outfile,title,colors):
	counts=df.groupby(["year","mlst_tar"]).size().unstack(fill_value=0)
	columns=counts.sum().sort_values(ascending=False).index.tolist()
	if "Other" in columns:
		columns.remove("Other")
		columns.append("Other")
	counts=counts.reindex(columns=columns,fill_value=0)
	percentages=counts.div(counts.sum(axis=1),axis=0)*100
	plot_colors=[]
	color_index=0
	for column in columns:
		if column=="Other":
			plot_colors.append("#808080")
		else:
			plot_colors.append(colors[color_index%len(colors)])
			color_index+=1
	fig,ax=plt.subplots(figsize=(11,7))
	ax.stackplot(
		percentages.index,
		*[percentages[column] for column in columns],
		labels=columns,
		colors=plot_colors
	)
	ax.set_ylabel("Percentage of genomes")
	ax.set_xlabel("Year")
	ax.set_ylim(0,100)
	ax.set_title(title)
	ax.set_xticks(sorted(percentages.index)[::2])
	ax.legend(title="MLST + tar",bbox_to_anchor=(1.02,1),loc="upper left",fontsize=8)
	ax.set_aspect("auto")
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close(fig)





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
			["year","mlst_tar","count"],
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
			["year","mlst_tar","count"],
		)
		return
	df,keep=prepare_data(
		df,
		smk.params.mlst_column,
		smk.params.date_column,
		float(smk.params.threshold),
	)
	if df.empty:
		if geography=="world":
			raise ValueError("No valid MLST/tar temporal data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No valid MLST/tar temporal data for {geography}: {geography_value}",
			["year","mlst_tar","count"],
		)
		return
	table=df.groupby(["year","mlst_tar"]).size().reset_index(name="count")
	table.to_csv(smk.output.table,sep="\t",index=False)
	plot_over_time(
		df,
		smk.output.plot,
		f"MLST + tar combination over time - {geography}"
		+("" if geography=="world" else f": {geography_value}")
		+f" ({smk.wildcards.dataset}, {smk.wildcards.interval})",
		smk.params.mlst_tar_colors
	)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore