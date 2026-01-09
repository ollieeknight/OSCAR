# from https://github.com/10XGenomics/cellranger/issues/133#issuecomment-989119662

#!/usr/bin/env python3

import argparse
import csv
from natsort import natsorted
import pandas as pd

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument(
  '-i',
  dest = 'input',
  required = True,
  help = 'Input gtf'
)
parser.add_argument(
  '-o',
  dest = 'output',
  required = True,
  help = 'Output sorted gtf'
)
args = parser.parse_args()

# Define the columns and their data types for the GTF file
gtf_columns = {
  'chromosome': 'str',
  'source': 'str',
  'feature': 'str',
  'start': 'uint64',
  'end': 'uint64',
  'score': 'str',
  'strand': 'str',
  'frame': 'str',
  'attribute': 'str'
}

# Read the GTF file into a pandas DataFrame
gtf = pd.read_csv(
  args.input,
  sep = '\t',
  comment = '#',
  names = gtf_columns.keys(),
  dtype = gtf_columns
)

# Extract gene_id and transcript_id from the attribute column
gtf['gene'] = gtf['attribute'].str.extract(r'gene_id "([^"]*)"')
gtf['transcript'] = gtf['attribute'].str.extract(r'transcript_id "([^"]*)"').fillna('')

# Get the start and end positions for each gene
genes = gtf[gtf['feature'] == 'gene'][
  ['chromosome', 'strand', 'gene', 'start', 'end']
].rename(
  columns = {
    'start': 'gene_start',
    'end': 'gene_end'
  }
).set_index(
  ['chromosome', 'strand', 'gene']
)

# Get the start and end positions for each transcript
transcripts = gtf[gtf['feature'] == 'transcript'][
  ['chromosome', 'strand', 'transcript', 'start', 'end']
].rename(
  columns = {
    'start': 'transcript_start',
    'end': 'transcript_end',
  }
).set_index(
  ['chromosome', 'strand', 'transcript']
)

# Add gene and transcript start and end positions to each row
gtf = gtf.set_index(
  ['chromosome', 'strand', 'gene']
).merge(
  genes,
  how = 'left',
  on = ['chromosome', 'strand', 'gene']
).reset_index().set_index(
  ['chromosome', 'strand', 'transcript']
).merge(
  transcripts,
  how = 'left',
  on = ['chromosome', 'strand', 'transcript']
).reset_index()

# Sort the rows based on multiple columns
gtf = gtf.sort_values(
  by = [
    'chromosome',
    'strand',
    'gene_start',
    'gene_end',
    'gene',
    'transcript_start',
    'transcript_end',
    'transcript',
    'feature',
    'start',
    'end'
  ],
  key = lambda x: (
    [0 if i == 'gene' else 1 if i == 'transcript' else 2 for i in x]
    if x.name == 'feature'
    else natsorted(x)
  )
)

# Write the sorted GTF to the output file
gtf.to_csv(
  args.output,
  sep = '\t',
  columns = list(gtf_columns.keys()),
  header = False,
  index = False,
  quoting = csv.QUOTE_MINIMAL,
  float_format = '%.10g'
)