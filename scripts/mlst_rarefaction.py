from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
	snakemake=None

import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline # type: ignore
from matplotlib.ticker import MaxNLocator # type: ignore
from utils import load_data,load_types


def get_sampling_sizes(total,step):
	sizes=np.arange(step,total+1,step)
	if len(sizes)==0 or sizes[-1]!=total:
		sizes=np.append(sizes,total)
	return sizes


def rarefy(values,sizes,repeats,seed):
	rng=np.random.default_rng(seed)
	results=[]
	for size in sizes:
		sampled_values=[]
		for _ in range(repeats):
			sampled=rng.choice(
				values,
				size=size,
				replace=False,
			)
			sampled_values.append(len(np.unique(sampled)))
		results.append({
			"genomes":size,
			"mean_STs":np.mean(sampled_values),
			"sd_STs":np.std(sampled_values),
		})
	return pd.DataFrame(results)


def run_rarefaction(values,step,repeats,seed):
	return rarefy(
		values,
		get_sampling_sizes(len(values),step),
		repeats,
		seed,
	)


def smooth_curve(result):
	x=result["genomes"].to_numpy()
	y=result["mean_STs"].to_numpy()
	if len(x)<4:
		return x,y
	x_smooth=np.linspace(x.min(),x.max(),500)
	y_smooth=make_interp_spline(
		x,
		y,
		k=3,
	)(x_smooth)
	return x_smooth,y_smooth


def grouped_rarefaction(
	df,
	group_column,
	mlst_column,
	step,
	repeats,
	seed,
):
	results={}
	for group,sub in df.groupby(group_column):
		values=load_types(sub,mlst_column)
		if len(values)<step:
			print(
				"SKIP",
				group,
				"because",
				len(values),
				"<",
				step,
			)
			continue
		results[str(group)]=run_rarefaction(
			values,
			step,
			repeats,
			seed,
		)
	return results


def plot_rarefaction(
	all_result,
	group_results,
	geography,
	outfile,
	title,
	geography_color_map,
):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(10,6))
	if all_result.empty:
		ax.axis("off")
		ax.text(0.5,0.5,"No data",ha="center",va="center")
		fig.savefig(outfile,bbox_inches="tight")
		plt.close(fig)
		return
	if geography=="world":
		x,y=smooth_curve(all_result)
		ax.plot(
			x,
			y,
			label="world",
			linewidth=2,
			color="#808080",
		)
	else:
		for name,result in group_results.items():
			x,y=smooth_curve(result)
			color=geography_color_map.get(str(name),"#808080")
			ax.plot(
				x,
				y,
				label=name,
				color=color,
			)
	ax.set_title(title)
	ax.set_xlabel("# genomes")
	ax.set_ylabel("Observed MLST types")
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	ax.legend(
		fontsize=8,
		bbox_to_anchor=(1.05,1),
		loc="upper left",
	)
	ax.grid(True,linestyle="--",linewidth=0.5,alpha=0.4)
	plt.tight_layout()
	fig.savefig(outfile,dpi=300,bbox_inches="tight")
	plt.close(fig)


def run_analysis(smk):
	df=load_data(smk.input.filtered)
	mlst_column=smk.params.mlst_column
	geography_column=smk.params.geography_column
	step=int(smk.params.step)
	repeats=int(smk.params.repeats)
	seed=int(smk.params.seed)
	all_values=load_types(df,mlst_column)
	if len(all_values)<step:
		raise ValueError(
			f"Not enough genomes for world: "
			f"{len(all_values)} < {step}"
		)
	all_result=run_rarefaction(
		all_values,
		step,
		repeats,
		seed,
	)
	if geography_column is None:
		group_results={}
	else:
		group_results=grouped_rarefaction(
			df,
			geography_column,
			mlst_column,
			step,
			repeats,
			seed,
		)
	all_result.to_csv(
		smk.output.table,
		sep="\t",
		index=False,
	)
	plot_rarefaction(
		all_result,
		group_results,
		smk.params.geography,
		smk.output.plot,
		f"MLST rarefaction - {smk.params.geography} "
		f"({smk.wildcards.dataset}, {smk.wildcards.interval})",
		smk.params.geography_color_map,
	)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore