process GenMasterdata {
	label 'process_low'
	publishDir params.out_dir, mode: 'copy'

	input:
	path known_list
	path ref_gene
	path ref_fasta
	val flexi_searchlen
	val fuscia_up
	val fuscia_down
	path gen_masterdata_py

	output:
	path "masterdata.csv"	

	script:
	"""
	python $gen_masterdata_py \
		$known_list \
		$ref_gene \
		$ref_fasta \
		$flexi_searchlen . $fuscia_up $fuscia_down
	"""
}
