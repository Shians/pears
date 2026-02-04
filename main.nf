// setting up scrips and software
include { GenMasterdata } from './subworkflows/gen_masterdata.nf'
include { RunSTARSolo } from './subworkflows/run_STARsolo.nf'
include { runFuscia } from './subworkflows/fuscia.nf'
include { runFlexiplex } from './subworkflows/flexiplex.nf'
include { runArriba } from './subworkflows/arriba.nf'
include { formatFuscia } from './subworkflows/formatting.nf'
include { getBarcodesArriba } from './subworkflows/arriba.nf'
include {formatFlexiplex as formatFlexiplex1} from './subworkflows/formatting.nf'
include {formatFlexiplex as formatFlexiplex2} from './subworkflows/formatting.nf'

workflow {
	masterdata_ch = GenMasterdata(
		file(params.known_list),
		file(params.ref_gene),
		file(params.ref_fasta),
		params.flexi_searchlen,
		params.fuscia_up,
		params.fuscia_down,
		file("$projectDir/subworkflows/gen_masterdata.py")
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

	STARsolo_result = RunSTARSolo(
		Channel.fromPath(params.read1).collect(),
		Channel.fromPath(params.read2).collect(),
		file(params.genome_index),
		file(params.barcode_whitelist),
		params.umi_len
	)

	Fuscia_output_ch   = runFuscia(mapped_ch, STARsolo_result[0], STARsolo_result[1])
	Flexiplex_output_ch = runFlexiplex(mapped_ch)
	Arriba_output_ch = runArriba(STARsolo_result[0])
	ArribaBC_output_ch  = getBarcodesArriba(mapped_ch, Arriba_output_ch)

	// collapse each into a single emission
	Fuscia_collected   = Fuscia_output_ch | collect
	Flexiplex_collected = Flexiplex_output_ch | collect
	ArribaBC_collected  = ArribaBC_output_ch | collect

	// formatting
	formatFuscia(Fuscia_collected, "master_fuscia.csv")
	formatFlexiplex1(Flexiplex_collected, "flexiplex_out", "master_flexiplex.csv")
	formatFlexiplex2(ArribaBC_collected, "arriba_out", "master_arriba.csv")
}
