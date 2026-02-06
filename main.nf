nextflow.enable.dsl=2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { downloadReferences } from './modules/download_references.nf'
include { genFusionTargets } from './modules/gen_fusion_targets.nf'
include { prepareIncludeList } from './modules/prepare_include_list.nf'
include { calculateReadLength } from './modules/calculate_read_length.nf'
include { buildSTARIndex; runSTARSolo } from './modules/star_solo.nf'
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

	// Use pre-built references if all three are provided, otherwise download
	if (params.ref_fasta && params.ref_gtf && params.star_genome_index) {
		log.info "Using pre-built references: skipping download and index building"
		ref_fasta = channel.value(file(params.ref_fasta))
		ref_gtf = channel.value(file(params.ref_gtf))
		star_index = channel.value(file(params.star_genome_index))
	} else {
		// Download reference genome and annotation
		references = downloadReferences(params.genome_version)
		ref_fasta = references.fasta
		ref_gtf = references.gtf
		star_index = null  // Will be built below
	}

	// Prepare barcode include list (decompress if gzipped)
 	include_list  = prepareIncludeList(file(params.barcode_include_list))

	fusion_targets = genFusionTargets(
		file(params.known_fusions_list),
		ref_gtf,
		ref_fasta,
		params.flexiplex_searchlen,
		params.fuscia_up,
		params.fuscia_down
	)

	fusion_target_rows = fusion_targets \
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

	// Build STAR index if not already provided via pre-built references
	if (!params.star_genome_index) {
		// Calculate R2 read length for STAR index generation
		r2_files = channel.fromPath(params.fastq_r2).collect()
		read_length = calculateReadLength(r2_files)

		star_index = buildSTARIndex(
			ref_fasta,
			ref_gtf,
			read_length
		)
	}

	STARsolo_result = runSTARSolo(
		channel.fromPath(params.fastq_r1).collect(),
		channel.fromPath(params.fastq_r2).collect(),
		star_index,
		include_list,
		params.umi_len
	)

	Fuscia_output   = runFuscia(fusion_target_rows, STARsolo_result.bam, STARsolo_result.bam_index, params.fuscia_mapqual)
	Flexiplex_output = runFlexiplex(
		fusion_target_rows,
		include_list,
		channel.fromPath(params.fastq_r1).collect(),
		channel.fromPath(params.fastq_r2).collect(),
		params.flexiplex_demultiplex_options
	)
	Arriba_output = runArriba(STARsolo_result.bam, ref_fasta, ref_gtf)
	ArribaBC_output  = getBarcodesArriba(
		fusion_target_rows,
		Arriba_output,
		include_list,
		channel.fromPath(params.fastq_r1).collect(),
		params.flexiplex_demultiplex_options
	)

	// collapse each into a single emission
	Fuscia_collected = Fuscia_output | collect
	Flexiplex_collected = Flexiplex_output | collect
	ArribaBC_collected = ArribaBC_output | collect

	// formatting
	formatFuscia(Fuscia_collected, "master_fuscia.csv")
	formatFlexiplex1(Flexiplex_collected, "flexiplex_out", "master_flexiplex.csv", params.out_dir)
	formatFlexiplex2(ArribaBC_collected, "arriba_out", "master_arriba.csv", params.out_dir)
}
