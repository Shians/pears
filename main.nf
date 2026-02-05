include { genMasterdata } from './modules/gen_masterdata.nf'
include { prepareIncludeList } from './modules/prepare_include_list.nf'
include { runSTARSolo } from './modules/run_STARsolo.nf'
include { runFuscia } from './modules/fuscia.nf'
include { runFlexiplex } from './modules/flexiplex.nf'
include { runArriba } from './modules/arriba.nf'
include { formatFuscia } from './modules/formatting.nf'
include { getBarcodesArriba } from './modules/arriba.nf'
include { formatFlexiplex as formatFlexiplex1 } from './modules/formatting.nf'
include { formatFlexiplex as formatFlexiplex2 } from './modules/formatting.nf'

workflow {
	// Prepare barcode include list (decompress if gzipped)
	include_list_ch = prepareIncludeList(file(params.barcode_include_list))

	masterdata_ch = genMasterdata(
		file(params.known_fusions_list),
		file(params.ref_gene),
		file(params.ref_fasta),
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

	STARsolo_result = runSTARSolo(
		channel.fromPath(params.fastq_r1).collect(),
		channel.fromPath(params.fastq_r2).collect(),
		file(params.star_genome_index),
		include_list_ch,
		params.umi_len
	)

	Fuscia_output_ch   = runFuscia(mapped_ch, STARsolo_result[0], STARsolo_result[1])
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
