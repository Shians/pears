process formatFuscia {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path fuscia_files
	val output_file

	output:
	path "${output_file}"

	script:
	"""
	#!/usr/bin/env python3

	import pandas as pd
	import os
	import glob

	def add_fusion_name(file):
		r = pd.read_table(file)
		if r.empty == False:
			r['fusion'] = os.path.basename(file).split("_")[0]
			r = r[r['cell_barcode'] != '-']
			return r

	df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'fusion'])
	for file in glob.glob('*'):
		if os.path.isfile(file):
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
	path('*')
	val dir
	val output_file

	output:
	path "${output_file}"

	script:
	"""
	#!/usr/bin/env python3

	import pandas as pd
	import os

	def add_fusion_name(file):
		r = pd.read_table(file)
		if r.empty == False:
			df_temp = pd.DataFrame().assign(cell_barcode = r['CellBarcode'], molecular_barcode = r['UMI'])
			df_temp['fusion'] =  os.path.basename(file).split("_")[1]
			return df_temp

	in_dir = f'$params.out_dir/$dir/'
	df = pd.DataFrame(columns = ['cell_barcode', 'molecular_barcode', 'fusion'])
	for file in os.listdir(in_dir):
		if os.path.basename(file)[0:8] == 'barcodes':
			r = add_fusion_name(f'{in_dir}/{file}')
			df = pd.concat([df, r], axis = 0, ignore_index = True)
	# Remove non-unique rows
	df = df.drop_duplicates()
	df.to_csv(f'$output_file', index=False)

	"""

}


