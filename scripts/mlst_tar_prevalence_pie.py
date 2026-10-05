from pathlib import Path
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import pandas as pd
from utils import load_data


TAR_CATEGORIES=["none","tarM","tarM+tarP","tarM+tarP+tarS","tarM+tarS","tarP","tarP+tarS","tarS"]


def write_empty(outfile,table,title,columns):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(10,8))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)
	pd.DataFrame(columns=columns).to_csv(table,sep="\t",index=False)


def prepare_mlst(df,mlst_column,threshold):
	df=df.dropna(subset=[mlst_column]).copy()
	counts=df[mlst_column].value_counts()
	cum=counts.cumsum()/counts.sum()
	keep=list(cum[cum<=threshold].index)
	if not keep:
		keep=[counts.index[0]]
	df["mlst_group"]=df[mlst_column].where(df[mlst_column].isin(keep),"Other")
	return df


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


def prepare_table(df):
	counts=df.groupby(["mlst_group","tar_combination"]).size().reset_index(name="count")
	mlst_totals=counts.groupby("mlst_group")["count"].sum()
	total=counts["count"].sum()
	counts["mlst_percentage"]=counts["count"]/counts["mlst_group"].map(mlst_totals)*100
	counts["percentage"]=counts["count"]/total*100
	return counts


def plot_pie(df,outfile,title,st_color_map,tar_colors):
	mlst_counts=df.groupby("mlst_group")["count"].sum()
	mlst_order=mlst_counts.sort_values(ascending=False).index.tolist()
	if "Other" in mlst_order:
		mlst_order.remove("Other")
		mlst_order.append("Other")
	tar_order=[category for category in TAR_CATEGORIES if category in df["tar_combination"].values]
	mlst_colors={}
	for mlst in mlst_order:
		if mlst=="Other":
			mlst_colors[mlst]="#808080"
		else:
			mlst_colors[mlst]=st_color_map.get(str(mlst),"#808080")
	inner_sizes=[mlst_counts[mlst] for mlst in mlst_order]
	inner_colors=[mlst_colors[mlst] for mlst in mlst_order]
	outer_sizes=[]
	outer_colors=[]
	for mlst in mlst_order:
		sub=df[df["mlst_group"]==mlst]
		for tar in tar_order:
			value=sub.loc[sub["tar_combination"]==tar,"count"].sum()
			if value>0:
				outer_sizes.append(value)
				outer_colors.append(tar_colors[tar])
	fig,ax=plt.subplots(figsize=(10,8))
	ax.pie(inner_sizes,radius=0.9,colors=inner_colors,startangle=90,wedgeprops={"edgecolor":"black"})
	ax.pie(outer_sizes,radius=0.9,colors=outer_colors,startangle=90,wedgeprops={"width":0.16,"edgecolor":"white"})
	mlst_percentages=mlst_counts/mlst_counts.sum()*100
	mlst_handles=[
		plt.Line2D([0],[0],marker="o",color="none",markerfacecolor=mlst_colors[mlst],markersize=8,label=f"{mlst} (n={mlst_counts[mlst]}, {mlst_percentages[mlst]:.1f}%)")
		for mlst in mlst_order
	]
	tar_counts=df.groupby("tar_combination")["count"].sum()
	tar_percentages=tar_counts/tar_counts.sum()*100
	tar_handles=[
		plt.Line2D([0],[0],marker="o",color="none",markerfacecolor=tar_colors[tar],markersize=8,label=f"{tar} (n={tar_counts[tar]}, {tar_percentages[tar]:.1f}%)")
		for tar in tar_order
	]
	legend1=ax.legend(handles=mlst_handles,title="MLST",bbox_to_anchor=(1.02,1),loc="upper left",fontsize=8)
	ax.text(0,-0.3,"MLST",ha="center",va="center",fontsize=11)
	ax.text(0,0.95,"tar",ha="center",va="center",fontsize=11)
	ax.add_artist(legend1)
	ax.legend(handles=tar_handles,title="Tar combination",bbox_to_anchor=(1.02,0.35),loc="upper left",fontsize=8)
	ax.set_title(title)
	ax.set(aspect="equal")
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	plt.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close()


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
			["mlst_group","tar_combination","count","mlst_percentage","percentage"],
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
			["mlst_group","tar_combination","count","mlst_percentage","percentage"],
		)
		return
	df=prepare_mlst(df,smk.params.mlst_column,float(smk.params.threshold))
	if df.empty:
		if geography=="world":
			raise ValueError("No MLST data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No MLST data for {geography}: {geography_value}",
			["mlst_group","tar_combination","count","mlst_percentage","percentage"],
		)
		return
	table=prepare_table(df)
	if table.empty:
		if geography=="world":
			raise ValueError("No MLST/tar data for world")
		write_empty(
			smk.output.plot,
			smk.output.table,
			f"No MLST/tar data for {geography}: {geography_value}",
			["mlst_group","tar_combination","count","mlst_percentage","percentage"],
		)
		return
	table.to_csv(smk.output.table,sep="\t",index=False)
	title=f"MLST prevalence and tar combinations - {geography}"
	if geography!="world":
		title+=f": {geography_value}"
	title+=f" - genomes: {smk.wildcards.dataset} ({smk.wildcards.interval})"
	plot_pie(table,smk.output.plot,title,smk.params.st_color_map,smk.params.tar_colors)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore