#!/usr/bin/env nextflow

params.ascpbin = "/rds/homes/n/nicholsz/.aspera/cli/bin/ascp"
params.publish = "/cephfs/covid/bham/nicholsz/artifacts/elan2"
params.dhmanifest = "/cephfs/covid/software/sam/dh/20200421/manifest.txt"
//params.k2db = "/ramdisk/kraken2db"
params.k2db = "/data/kraken2db"

Channel
    .fromPath(params.manifest)
    .splitCsv(header:true, sep:'\t')
    .map { row-> tuple(row.ena_sample_name, row.central_sample_id, row.sample_center_name, row.collection_date, row.received_date, row.adm0, row.adm1, row.run_name, row.published_name, file(row.climb_fn), row.run_center_name, row.library_strategy, row.library_source, row.library_selection, row.instrument_make, row.instrument_model, row.virus_identifier, row.min_ct_value, row.max_ct_value, row.library_primers, row.library_protocol, row.library_seq_kit, row.library_seq_protocol) }
    .set { manifest_ch }

process extract_bam_reads {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from manifest_ch

    output:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file("${bam}.fasta"), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol into bamfa_manifest_ch

    """
    samtools view ${bam} | awk '{print ">"\$1"\\n"\$10}' > ${bam}.fasta
    """
}

process kraken_bam_reads {
    tag { bam }
    conda "environments/kraken2.yaml"

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from bamfa_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/k2", pattern: "*k2*", mode: "copy", overwrite: true
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file("${bam_fasta}.k2o.9606.ls"), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol into k2_manifest_ch
    file "${bam_fasta}.k2o"
    file "${bam_fasta}.k2r"

    cpus 4
    """
    kraken2 --db ${params.k2db} --threads ${task.cpus} --output ${bam_fasta}.k2o --report ${bam_fasta}.k2r ${bam_fasta} && awk '\$3 == 9606 {print \$2}' ${bam_fasta}.k2o > ${bam_fasta}.k2o.9606.ls
    """
}

process dehumanise_bam {
    tag { bam }
    conda "environments/dehumanizer.yaml"
    label 'bear'

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from k2_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/dh", pattern: "${coguk_id}.${run_name}.dh", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/staging/alignment-clean/", pattern: "${coguk_id}.${run_name}.climb.public.bam", mode: "copy", overwrite: true
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), file("${coguk_id}.${run_name}.climb.public.bam"), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol into ascp_manifest_ch
    file "${coguk_id}.${run_name}.climb.public.bam"
    file "${coguk_id}.${run_name}.dh" into dh_report_ch

    errorStrategy 'retry' 
    maxRetries 3
    memory { (12 + (2 * task.attempt))+"GB" }

    script:
    if ( run_platform == "ILLUMINA" )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset sr --bam -o ${coguk_id}.${run_name}.climb.public.bam --trash-minalen 25 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls}
        """
    else if( run_platform == 'OXFORD_NANOPORE' )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset map-ont --bam -o ${coguk_id}.${run_name}.climb.public.bam --trash-minalen 10 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls}
        """
    else
        error "Invalid alignment mode for technology ${run_platform}"
}

process ascp_bam {
    tag { bam }
    cpus 6 //# massively over-request local cores to prevent sending too much to API at once
    errorStrategy { sleep(Math.pow(2, task.attempt) * 300 as long); return 'retry' }
    maxRetries 5

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), file(public_bam), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from ascp_manifest_ch
    output:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), file(public_bam), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol into publish_manifest_ch

    """
    export ASPERA_SCP_PASS=\$WEBIN_PASS
    ${params.ascpbin} -T --policy high -L- ${public_bam} \$WEBIN_USER@webin.ebi.ac.uk:.
    """
}

process publish_bam {
    tag { bam }
    conda "environments/pyena.yaml"

    cpus 6 //# massively over-request local cores to prevent sending too much to API at once

    errorStrategy { sleep(Math.pow(2, task.attempt) * 300 as long); return 'retry' }
    maxRetries 1

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), file(public_bam), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from publish_manifest_ch

    output:
    file "${public_bam}.txt" into dh_accession_report_ch
    file "${public_bam}.txt" into dh_ocarina_report_ch

    script:
    """
    pyena --study-accession ${params.study} --my-data-is-ready --no-ftp \
          --sample-name ${ena_sample_name} \
          --sample-center-name '${sample_center}' \
          --sample-taxon '2697049' \
          --sample-attr 'collector name' 'not provided' \
          --sample-attr 'collecting institution' '${sample_center}' \
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
          --run-center-name '${run_center}' \
          --run-instrument '${run_instrument}' \
          --run-lib-protocol '${exp_seq_kit}|${exp_seq_protocol}' \
          --run-lib-source ${l_source} \
          --run-lib-selection ${l_selection} \
          --run-lib-strategy ${l_strategy} > ${public_bam}.txt
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

    cpus 6 //# massively over-request local cores to prevent sending too much to API at once

    script:
    """
    ocarina --env put publish --publish-group '${ena_run_name}' --service 'ENA-SAMPLE' --accession ${sample_acc} --public
    ocarina --env put publish --publish-group '${ena_run_name}' --service 'ENA-RUN' --accession ${run_acc} --public
    """
}

dh_report_ch
    .collectFile(name: "dehumanised.qc", storeDir: "${params.publish}/staging/dh/${params.datestamp}/", keepHeader: true)

dh_accession_report_ch
    .collectFile(name: "dehumanised.accessions.ls", storeDir: "${params.publish}/staging/dh/${params.datestamp}/")

