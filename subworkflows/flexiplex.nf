//FIX -d option set in config file.


//path "${fusion_genes}_${task.index}_reads.fastq"
//path barcode_file

process runFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy'
	
	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)

	output:
	path "barcodes_${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_reads_barcodes.txt"

	script:
	fusion_name="${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"
	"""
	# Define the fusion name and flexiplex path

	# Run flexiplex with the specified parameters
	paste <(gunzip -c ${params.fastq_r1}) <(gunzip -c ${params.fastq_r2}) | \
	sed "/^[@+]/! s/^/START/g" | sed "/^[@+]/! s/	//g" | \
	${projectDir}/modules/flexiplex/flexiplex -p $task.cpus -n ${fusion_name} \
		-x ${sequence1}${sequence2} -d grep -f 1 > ${fusion_name}_reads.fastq

	${projectDir}/modules/flexiplex/flexiplex -x START \
		${params.flexiplex_demultiplex_options} \
		-k ${params.barcode_whitelist} -n barcodes_${fusion_name} ${fusion_name}_reads.fastq
    """
}

