from pathlib import Path

import matplotlib # type: ignore
matplotlib.use("Agg")

import matplotlib.pyplot as plt # type: ignore
import pandas as pd
import geopandas as gpd # type: ignore
import numpy as np

from matplotlib.patches import Wedge # type: ignore

from utils import load_data, get_mrsa


def prepare_mrsa(df, region_column, mrsa_column, mrsa_value):

	df=df.dropna(subset=[region_column]).copy()

	df["MRSA"]=get_mrsa(
		df,
		mrsa_column,
		mrsa_value
	)

	out=(
		df.groupby([region_column,"MRSA"])
		.size()
		.reset_index(name="count")
	)

	pivot=(
		out
		.pivot(
			index=region_column,
			columns="MRSA",
			values="count"
		)
		.fillna(0)
	)

	pivot["total"]=pivot.sum(axis=1)

	pivot["MRSA"]=pivot.get(True,0)
	pivot["MSSA"]=pivot.get(False,0)

	pivot=pivot.reset_index()

	return pivot


def load_europe_map(df, region_column):

	world=gpd.read_file(
		gpd.datasets.get_path("naturalearth_lowres")
	)

	europe=world[
		world["continent"]=="Europe"
	].copy()

	# rename to match your metadata
	europe=europe.merge(
		df[["country",region_column]]
		.drop_duplicates(),
		left_on="name",
		right_on="country",
		how="left"
	)

	europe=europe.dropna(
		subset=[region_column]
	)

	# create UNRegion polygons
	regions=europe.dissolve(
		by=region_column
	)

	return europe, regions


def draw_pie(ax, x, y, mrsa, mssa, radius):

	total=mrsa+mssa

	if total==0:
		return

	start=0

	for value,color in [
		(mrsa,"red"),
		(mssa,"steelblue")
	]:

		angle=value/total*360

		wedge=Wedge(
			(x,y),
			radius,
			start,
			start+angle,
			facecolor=color,
			edgecolor="white",
			linewidth=0.5,
			alpha=0.85
		)

		ax.add_patch(wedge)

		start+=angle


def plot_map(
	df,
	europe,
	regions,
	region_column,
	outfile
):

	fig,ax=plt.subplots(
		figsize=(8,8)
	)


	# countries
	europe.plot(
		ax=ax,
		color="white",
		edgecolor="lightgrey",
		linewidth=0.3
	)


	# UNRegion borders
	regions.boundary.plot(
		ax=ax,
		color="black",
		linewidth=1
	)


	# region coordinates
	centroids=regions.geometry.centroid

	for _,row in df.iterrows():

		if row[region_column] not in centroids:
			continue

		point=centroids[row[region_column]]

		# scale pie size
		radius=(
			np.log2(row["total"]+1)-2
		)/100

		draw_pie(
			ax,
			point.x,
			point.y,
			row["MRSA"],
			row["MSSA"],
			radius
		)


	ax.set_xlim(
		europe.total_bounds[0],
		europe.total_bounds[2]
	)

	ax.set_ylim(
		europe.total_bounds[1],
		europe.total_bounds[3]
	)

	ax.axis("off")


	# legend
	ax.scatter( [], [], color="red", label="MRSA" )

	ax.scatter( [], [], color="steelblue", label="MSSA" )

	ax.legend( frameon=False, loc="lower left" )


	Path(outfile).parent.mkdir( parents=True, exist_ok=True )

	plt.savefig( outfile, dpi=300, bbox_inches="tight" )

	plt.close()


def run_analysis(smk):

	df=load_data(
		smk.input.filtered
	)

	df_summary=prepare_mrsa(
		df,
		smk.params.region_column,
		smk.params.mrsa_column,
		smk.params.mrsa_value
	)


	europe,regions=load_europe_map(
		df,
		smk.params.region_column
	)


	plot_map(
		df_summary,
		europe,
		regions,
		smk.params.region_column,
		smk.output.plot
	)


	df_summary.to_csv(
		smk.output.table,
		sep="\t",
		index=False
	)


if __name__=="__main__":
	run_analysis(snakemake) # type: ignore