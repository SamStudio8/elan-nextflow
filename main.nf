#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {play_ocarina} from "./modules/ocarina.nf"
include {inbound} from "./workflows/elan.nf"

if (!params.mode) error "A workflow must be chosen: --mode {inbound,ocarina,ena_bam}"
if (!params.ocarina_profile) error "An Ocarina profile must be chosen: --ocarina_profile profile"
if( !System.getenv("OCARINA_CONF_FILE") ) error '$OCARINA_CONF_FILE unset'

workflow {

    if (params.mode == "inbound"){
        if( !params.datestamp ) error "Missing `datestamp` param: YYYYMMDD datestamp to identify today's run"
        if( !params.uploads ) error "Missing `uploads` param: path to glob CLIMB-COVID uploads"
        if( !params.uploads_usern ) error "Missing `uploads_usern` param: index of username directory in uploads glob, root counts as zero [int]"
        if( !params.publish ) error "Missing `publish` param: path to CLIMB-COVID staged artifacts root"
        if( !params.minlen ) error "Missing `min_len` param: minimum genome size required to pass the save_uploads step [int]"
        if( !params.artifacts_root ) error "Missing `artifacts_root` param: path to new CLIMB-COVID published artifacts root"

        if( !System.getenv("ELAN_SLACK_HOOK") ) error '$ELAN_SLACK_HOOK unset'

	inbound()
        
    } else if (params.mode == "ocarina"){
        if( !params.manifest ) error "Missing `manifest` param: today's ocarina manifest from elan [path]"
        
        Channel
            .fromPath(params.manifest)
            .splitCsv(header:['coguk_id', 'run_name', 'username', 'pipeuuid', 'dir', 'fasta', 'bam', 'qc', 'sourcesite', 'seqsite', 'platform'], sep:'\t')
            .map{ it << ["qc": [it["dir"], it["qc"]].join('/')] }
            .set{manifest_ch}

        play_ocarina(manifest_ch)

    } else if (params.mode == "ena_bam") {

    } else {
        error "Invalid mode selected"
    }
}
