#!/usr/bin/env nextflow

params.ascpbin = "/rds/homes/n/nicholsz/.aspera/cli/bin/ascp"
params.publish = "/cephfs/covid/bham/nicholsz/artifacts/elan2"

Channel
    .fromPath(params.manifest)
    .splitCsv(header:['ena_sample_name', 'coguk_id', 'sample_center', 'collection_date', 'received_date', 'adm0', 'adm1', 'run_name', 'ena_run_name', 'bam', 'run_center', 'l_strategy', 'l_source', 'l_selection', 'run_platform', 'run_instrument', 'bam_fasta', 'bam_hum_ls', 'public_bam', 'gisaid_id', 'min_ct', 'max_ct', 'exp_primers', 'exp_protocol', 'exp_seq_kit', 'exp_seq_protocol'], sep:'\t')
    .map { row-> tuple(row.ena_sample_name, row.coguk_id, row.sample_center, row.collection_date, row.received_date, row.adm0, row.adm1, row.run_name, row.ena_run_name, row.bam, row.run_center, row.l_strategy, row.l_source, row.l_selection, row.run_platform, row.run_instrument, row.bam_fasta, row.bam_hum_ls, [params.publish, 'staging', 'alignment-clean', row.public_bam].join('/'), row.gisaid_id, row.min_ct, row.max_ct, row.exp_primers, row.exp_protocol, row.exp_seq_kit, row.exp_seq_protocol) }
    .set { ascp_manifest_ch }

process ascp_bam {
    tag { bam }
    cpus 6 //# massively over-request local cores to prevent too many ascp jobs
    errorStrategy { sleep(Math.pow(2, task.attempt) * 300 as long); return 'retry' }
    maxRetries 5

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, bam, run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, bam_fasta, bam_hum_ls, public_bam, gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from ascp_manifest_ch
    output:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, bam, run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, bam_fasta, bam_hum_ls, public_bam, gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol into publish_manifest_ch

    """
    export ASPERA_SCP_PASS=\$WEBIN_PASS
    ${params.ascpbin} -T --policy high -L- ${public_bam} \$WEBIN_USER@webin.ebi.ac.uk:.
    """
}

process publish_bam {
    tag { bam }
    conda "environments/pyena.yaml"

    cpus 6 //# massively over-request local cores to prevent sending too much to API at once

    //errorStrategy { sleep(Math.pow(2, task.attempt) * 300 as long); return 'retry' }
    //maxRetries 1
    errorStrategy 'ignore'

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, bam, run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, bam_fasta, bam_hum_ls, public_bam, gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from publish_manifest_ch

    output:
    file "${coguk_id}.${run_name}.pyena.txt" into dh_accession_report_ch
    file "${coguk_id}.${run_name}.pyena.txt" into dh_ocarina_report_ch

    script:
    """
    pyena --study-accession ${params.study} --my-data-is-ready --no-ftp \
          --sample-name ${ena_sample_name} \
          --sample-center-name "${sample_center}" \
          --sample-taxon '2697049' \
          --sample-attr 'collector name' 'not provided' \
          --sample-attr 'collecting institution' "${sample_center}" \
          --sample-attr 'collection date' ${collection_date} \
          --sample-attr 'geographic location (country and/or sea)' 'United Kingdom' \
          --sample-attr 'geographic location (region and locality)' '${adm1}' \
          --sample-attr 'definition for seropositive sample' 'not provided' \
          --sample-attr 'serotype (required for a seropositive sample)' 'not provided' \
          --sample-attr 'host common name' 'not provided' \
          --sample-attr 'host health state' 'not provided' \
          --sample-attr 'host scientific name' 'Human' \
          --sample-attr 'host sex' 'not provided' \
          --sample-attr 'host subject id' 'not provided' \
          --sample-attr 'isolate' 'not provided' \
          --sample-attr 'receipt date' '${received_date}' \
          --sample-attr 'sample capture status' 'active surveillance in response to outbreak' \
          --sample-attr 'virus identifier' 'not provided' \
          --sample-attr 'ENA-CHECKLIST' 'ERC000033' \
          --sample-attr 'min_cycle_threshold' '${min_ct}' \
          --sample-attr 'max_cycle_threshold' '${max_ct}' \
          --experiment-attr 'artic_primer_version' '${exp_primers}' \
          --experiment-attr 'artic_protocol_version' '${exp_protocol}' \
          --run-name ${ena_run_name} \
          --run-file-path ${public_bam} \
          --run-center-name "${run_center}" \
          --run-instrument '${run_instrument}' \
          --run-lib-protocol '${exp_seq_kit}|${exp_seq_protocol}' \
          --run-lib-source ${l_source} \
          --run-lib-selection ${l_selection} \
          --run-lib-strategy ${l_strategy} > ${coguk_id}.${run_name}.pyena.txt
    """
}

dh_ocarina_report_ch
    .splitCsv(header:['success', 'real', 'ena_sample_name', 'ena_run_name', 'bam', 'study_acc', 'sample_acc', 'exp_acc', 'run_acc'], sep:' ')
    .map { row-> tuple(row.ena_run_name, row.sample_acc, row.run_acc) }
    .set { dh_ocarina_report_ch_split }

process tag_ocarina {
    tag { bam }
    label 'ocarina'
    conda "environments/ocarina.yaml"

    input:
    tuple ena_run_name, sample_acc, run_acc from dh_ocarina_report_ch_split

    errorStrategy { sleep(Math.pow(2, task.attempt) * 300 as long); return 'retry' }
    maxRetries 3

    cpus 6 //# massively over-request local cores to prevent sending too much to API at once

    script:
    """
    ocarina --oauth --env put publish --publish-group '${ena_run_name}' --service 'ENA-SAMPLE' --accession ${sample_acc} --public --submitted
    ocarina --oauth --env put publish --publish-group '${ena_run_name}' --service 'ENA-RUN' --accession ${run_acc} --public --submitted
    """
}

dh_accession_report_ch
    .collectFile(name: "dehumanised.accessions.ls", storeDir: "${params.publish}/staging/dh/${params.datestamp}/")

