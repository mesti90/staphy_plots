import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

def format_mlst(df, column="mlst"):

	df = df.copy()

	def convert(x):
		if pd.isna(x):
			return x

		try:
			x = float(x)
			if x.is_integer():
				return f"ST{int(x)}"
		except:
			pass

		return x

	df[column] = df[column].apply(convert)

	return df


def load_data(input_file, mlst_colname=None):

	df = pd.read_csv(
		input_file,
		sep="\t"
	)

	if mlst_colname and mlst_colname in df.columns:
		df = format_mlst(df, mlst_colname)

	return df


def extract_year(df, date_column="Date"):

	if date_column not in df.columns:
		raise ValueError(
			f"Missing date column: {date_column}"
		)

	df = df.copy()

	date = df[date_column].astype(str)

	valid = date.str.match(
		r"^\d{4}(?:-\d{2}(?:-\d{2})?)?$"
	)

	df["year"] = (
		date.where(valid)
		.str[:4]
		.astype(float)
	)

	return df


def filter_year(df, minyear=None, maxyear=None, date_column="Date"):

	df = extract_year(
		df,
		date_column
	)

	if minyear is not None:
		df = df[df["year"] >= minyear]

	if maxyear is not None:
		df = df[df["year"] <= maxyear]

	return df.copy()


def load_types(df, column):

	df = df.dropna(
		subset=[column]
	)

	return df[column].astype(str).values


def get_mrsa(df, mrsa_column, mrsa_value):

	if mrsa_column not in df.columns:
		raise KeyError(
			f"Missing MRSA column: {mrsa_column}"
		)

	return (
		df[mrsa_column]
		.astype(str)
		== str(mrsa_value)
	)

def write_empty_plot(outfile,title,table=None,columns=None):
	Path(outfile).parent.mkdir(parents=True,exist_ok=True)
	fig,ax=plt.subplots(figsize=(10,6))
	ax.axis("off")
	ax.text(0.5,0.5,title,ha="center",va="center",fontsize=12)
	fig.savefig(outfile,bbox_inches="tight")
	plt.close(fig)
	if table is not None:
		import pandas as pd
		pd.DataFrame(columns=columns or []).to_csv(table,sep="\t",index=False)