#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {play_ocarina} from "./modules/ocarina.nf"

if (!params.mode) error "A workflow must be chosen: --mode {inbound,ocarina,ena_bam}"

workflow {
    if (params.mode == "inbound"){

    } else if (params.mode == "ocarina"){
        if( !params.manifest ) error "Missing `manifest` param: today's ocarina manifest from elan [path]"
        if( !System.getenv("MAJORA_DOMAIN") ) error '$MAJORA_DOMAIN unset, Majora credentials likely not loaded into environment' // just check for MAJORA_DOMAIN here
        
        Channel
            .fromPath(params.manifest)
            .splitCsv(header:['coguk_id', 'run_name', 'username', 'pipeuuid', 'dir', 'fasta', 'bam', 'qc', 'sourcesite', 'seqsite', 'platform'], sep:'\t')
            // .map { row-> tuple(row.coguk_id, row.run_name, row.username, row.pipeuuid, row.fasta, row.bam, [row.dir, row.qc].join('/'), row.sourcesite, row.seqsite, row.platform) }
            .set { manifest_ch }

        play_ocarina(
            manifest_ch.row.coguk_id, 
            manifest_ch.row.run_name,
            manifest_ch.row.username,
            manifest_ch.row.pipeuuid,
            manifest_ch.row.fasta,
            manifest_ch.row.bam,
            [manifest_ch.row.dir, manifest_ch.row.qc].join('/'),
            manifest_ch.row.sourcesite,
            manifest_ch.row.seqsite,
            manifest_ch.row.platform
        )
    } else if (params.mode == "ena_bam") {

    } else {
        error "Invalid mode selected"
    }
}
