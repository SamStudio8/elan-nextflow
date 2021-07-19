#!/usr/bin/env nextflow

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
    conda "../environments/samtools.yaml"
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
    conda "../environments/kraken2.yaml"

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from bamfa_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/k2", pattern: "*k2*", mode: "copy", overwrite: true
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file("${bam_fasta}.k2o.9606.ls"), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol into k2_manifest_ch
    file "${bam_fasta}.k2o.gz"
    file "${bam_fasta}.k2r"

    cpus 4
    """
    kraken2 --db ${params.k2db} --threads ${task.cpus} --output ${bam_fasta}.k2o --report ${bam_fasta}.k2r ${bam_fasta} && awk '\$3 == 9606 {print \$2}' ${bam_fasta}.k2o > ${bam_fasta}.k2o.9606.ls
    gzip ${bam_fasta}.k2o
    """
}

process dehumanise_bam {
    tag { bam }
    conda "../environments/dehumanizer.yaml"
    label 'bear'

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from k2_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/dh", pattern: "${coguk_id}.${run_name}.dh", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/staging/alignment-clean/", pattern: "${coguk_id}.${run_name}.climb.public.bam", mode: "copy", overwrite: true
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), file("${coguk_id}.${run_name}.climb.public.bam"), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol into ascp_manifest_ch
    file "${coguk_id}.${run_name}.climb.public.bam"
    file "${coguk_id}.${run_name}.dh" into dh_report_ch

    errorStrategy 'ignore'
    memory "24GB"

    script:
    if ( run_platform == "ILLUMINA" )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset sr --bam -o ${coguk_id}.${run_name}.climb.public.bam --trash-minalen 25 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls} --pg-date ${params.datestamp}
        """
    else if( run_platform == 'OXFORD_NANOPORE' )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset map-ont --bam -o ${coguk_id}.${run_name}.climb.public.bam --trash-minalen 10 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls} --pg-date ${params.datestamp}
        """
    else if( run_platform == 'ION_TORRENT' )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset sr --bam -o ${coguk_id}.${run_name}.climb.public.bam --trash-minalen 25 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls} --pg-date ${params.datestamp}
        """
    else
        error "Invalid alignment mode for technology ${run_platform}"
}

process ascp_ls {

    input:
    tuple ena_sample_name, coguk_id, sample_center, collection_date, received_date, adm0, adm1, run_name, ena_run_name, file(bam), run_center, l_strategy, l_source, l_selection, run_platform, run_instrument, file(bam_fasta), file(bam_hum_ls), file(public_bam), gisaid_id, min_ct, max_ct, exp_primers, exp_protocol, exp_seq_kit, exp_seq_protocol from ascp_manifest_ch

    output:
    file "${coguk_id}.${run_name}.ascp.ls" into ascp_list_ch

    """
    echo "${ena_sample_name}\t${coguk_id}\t${sample_center}\t${collection_date}\t${received_date}\t${adm0}\t${adm1}\t${run_name}\t${ena_run_name}\t${bam}\t${run_center}\t${l_strategy}\t${l_source}\t${l_selection}\t${run_platform}\t${run_instrument}\t${bam_fasta}\t${bam_hum_ls}\t${public_bam}\t${gisaid_id}\t${min_ct}\t${max_ct}\t${exp_primers}\t${exp_protocol}\t${exp_seq_kit}\t${exp_seq_protocol}" > ${coguk_id}.${run_name}.ascp.ls
    """
}

dh_report_ch
    .collectFile(name: "dehumanised.qc", storeDir: "${params.publish}/staging/dh/${params.datestamp}/", keepHeader: true)

ascp_list_ch
    .collectFile(name: "ascp.files.ls", storeDir: "${params.publish}/staging/dh/${params.datestamp}/")
