process formatFuscia {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(fuscia_files)
	val(output_file)

	output:
	path "${output_file}"

	script:
	def input_files = fuscia_files.collect { f -> f.name }.join(' ')
	"""
	#!/usr/bin/env python3

	import pandas as pd
	import os

	input_files = '${input_files}'.split()

	def add_fusion_name(file):
		r = pd.read_table(file)
		if r.empty == False:
			r['fusion'] = os.path.basename(file).split("_")[0]
			r = r[r['cell_barcode'] != '-']
			return r

	df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'fusion'])
	for file in input_files:
		r = add_fusion_name(file)
		if r is not None:
			df = pd.concat([df, r], axis=0, ignore_index=True)
	df['cell_barcode'] = df['cell_barcode'].str.replace('-1', "")

	# Remove the last two columns
	df = df.iloc[:, :-3]
	# Remove non-unique rows
	df = df.drop_duplicates()

	df.to_csv('${output_file}', index=False)
	"""
}

process formatFlexiplex {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(barcode_files)
	val(output_file)

	output:
	path "${output_file}"

	script:
	def input_files = barcode_files.collect { f -> f.name }.join(' ')
	"""
	#!/usr/bin/env python3

	import pandas as pd
	import os

	input_files = '${input_files}'.split()

	def add_fusion_name(file):
		r = pd.read_table(file)
		if r.empty == False:
			df_temp = pd.DataFrame().assign(cell_barcode = r['CellBarcode'], molecular_barcode = r['UMI'])
			df_temp['fusion'] = os.path.basename(file).split("_")[1]
			return df_temp

	df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'fusion'])
	for file in input_files:
		if os.path.basename(file).startswith('barcodes'):
			r = add_fusion_name(file)
			if r is not None:
				df = pd.concat([df, r], axis=0, ignore_index=True)
	# Remove non-unique rows
	df = df.drop_duplicates()
	df.to_csv('${output_file}', index=False)
	"""
}

process formatArriba {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(barcode_files)
	val(output_file)

	output:
	path "${output_file}"

	script:
	def input_files = barcode_files.collect { f -> f.name }.join(' ')
	"""
	#!/usr/bin/env python3

	import pandas as pd
	import os

	input_files = '${input_files}'.split()

	def add_fusion_name(file):
		r = pd.read_table(file)
		if r.empty == False:
			df_temp = pd.DataFrame().assign(cell_barcode = r['CellBarcode'], molecular_barcode = r['UMI'])
			df_temp['fusion'] = os.path.basename(file).split("_")[1]
			return df_temp

	df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'fusion'])
	for file in input_files:
		if os.path.basename(file).startswith('barcodes'):
			r = add_fusion_name(file)
			if r is not None:
				df = pd.concat([df, r], axis=0, ignore_index=True)
	# Remove non-unique rows
	df = df.drop_duplicates()
	df.to_csv('${output_file}', index=False)
	"""
}


