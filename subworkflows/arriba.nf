process runArriba {
       	publishDir "${params.out_dir}/arriba_out", mode: 'copy'
	time = '1d'
	memory = '150'

	input:
	path bam_file

	output:
	path "fusions.tsv"

	script:
	"""

	${projectDir}/modules/arriba/arriba \
	    -x $bam_file \
	    -o fusions.tsv \
	    -O fusions.discarded.tsv \
	    -a $params.ref_fasta \
	    -g $params.ref_gene \
	    -f blacklist
	"""
}


//path "${fusion_genes}_${task.index}.fastq"
//path barcodes_file

process getBarcodes_Arriba {
	echo = true
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

	input:
        tuple val (fusion_genes), val (chrom1), val (gene1), val (base1), val (sequence1), val (chrom2), val (gene2), val (base2), val (sequence2)
	path fusion_table
	params.possible_barcode_list

       	output:
	path "barcodes_${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_reads_barcodes.txt"

	shell:
	"""
	set +e

	fus=`echo !{fusion_genes} | sed 's/--/\t/g'` ;
	pos=`echo -e "!{chrom1}:!{base1}\t!{chrom2}:!{base2}"`
	fusion_name=`echo !{fusion_genes}_!{chrom1}_!{base1}_!{chrom2}_!{base2}`

	grep -e "\$fus" -e "\$pos" !{fusion_table} |\
	   cut -f30 |\
	   sed 's/,/ \\n/g' |\
	   sed 's/^/^@/g' |\
	   grep -f - <(gunzip -c !{params.read1}) -A3 --no-group-separator |\
           sed "/^[@+]/! s/^/START/g" > "\$fusion_name".fastq ;

	   ${projectDir}/modules/flexiplex/flexiplex -x START \
	   ${params.flexiplex_demultiplex_options} \
	   -k ${params.possible_barcode_list} -n barcodes_"\$fusion_name" \
	   "\$fusion_name".fastq

	"""

}

