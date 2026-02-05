include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { downloadReferences } from './modules/download_references.nf'
include { genMasterdata } from './modules/gen_masterdata.nf'
include { prepareIncludeList } from './modules/prepare_include_list.nf'
include { calculateReadLength } from './modules/calculate_read_length.nf'
include { buildSTARIndex; runSTARSolo } from './modules/run_STARsolo.nf'
include { runFuscia } from './modules/fuscia.nf'
include { runFlexiplex } from './modules/flexiplex.nf'
include { runArriba } from './modules/arriba.nf'
include { formatFuscia } from './modules/formatting.nf'
include { getBarcodesArriba } from './modules/arriba.nf'
include { formatFlexiplex as formatFlexiplex1 } from './modules/formatting.nf'
include { formatFlexiplex as formatFlexiplex2 } from './modules/formatting.nf'

workflow {
	// Validate parameters against schema
	validateParameters()
	log.info paramsSummaryLog(workflow)

	// Download reference genome and annotation
	references = downloadReferences(params.genome_version)
	ref_fasta = references.fasta
	ref_gtf = references.gtf

	// Prepare barcode include list (decompress if gzipped)
	include_list_ch = prepareIncludeList(file(params.barcode_include_list))

	masterdata_ch = genMasterdata(
		file(params.known_fusions_list),
		ref_gtf,
		ref_fasta,
		params.flexiplex_searchlen,
		params.fuscia_up,
		params.fuscia_down
	)

	mapped_ch = masterdata_ch \
		| splitCsv(header:true) \
		| map { row ->
			tuple(
				row.fusion_genes,
				row.chrom1,
				row.gene1,
				row.base1,
				row.sequence1,
				row.chrom2,
				row.gene2,
				row.base2,
				row.sequence2
			)
		}

	// Calculate R2 read length for STAR index generation
	r2_files_ch = channel.fromPath(params.fastq_r2).collect()
	read_length_ch = calculateReadLength(r2_files_ch)

	// Build STAR index if not provided, otherwise use existing
	if (params.star_genome_index) {
		star_index_ch = channel.value(file(params.star_genome_index))
	} else {
		star_index_ch = buildSTARIndex(
			ref_fasta,
			ref_gtf,
			read_length_ch
		)
	}

	STARsolo_result = runSTARSolo(
		channel.fromPath(params.fastq_r1).collect(),
		channel.fromPath(params.fastq_r2).collect(),
		star_index_ch,
		include_list_ch,
		params.umi_len
	)

	Fuscia_output_ch   = runFuscia(mapped_ch, STARsolo_result.out.bam, STARsolo_result.out.bam_index)
	Flexiplex_output_ch = runFlexiplex(mapped_ch, include_list_ch)
	Arriba_output_ch = runArriba(STARsolo_result[0])
	ArribaBC_output_ch  = getBarcodesArriba(mapped_ch, Arriba_output_ch, include_list_ch)

	// collapse each into a single emission
	Fuscia_collected = Fuscia_output_ch | collect
	Flexiplex_collected = Flexiplex_output_ch | collect
	ArribaBC_collected = ArribaBC_output_ch | collect

	// formatting
	formatFuscia(Fuscia_collected, "master_fuscia.csv")
	formatFlexiplex1(Flexiplex_collected, "flexiplex_out", "master_flexiplex.csv")
	formatFlexiplex2(ArribaBC_collected, "arriba_out", "master_arriba.csv")
}
