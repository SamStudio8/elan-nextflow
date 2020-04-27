#!/usr/bin/env nextflow

params.publish = "/cephfs/covid/bham/nicholsz/artifacts/elan2"
params.dhmanifest = "/cephfs/covid/software/sam/dh/20200421/manifest.txt"
//params.k2db = "/ramdisk/kraken2db"
params.k2db = "/data/kraken2db"

Channel
    .fromPath(params.manifest)
    .splitCsv(header:['central_sample_id', 'instrument_make', 'run_name'], sep:'\t')
    .map { row-> tuple(row.instrument_make, row.central_sample_id, row.run_name, file([params.publish, 'staging', 'alignment', "${row.central_sample_id}.${row.run_name}.climb.bam"].join('/'))) }
    .set { manifest_ch }

process extract_bam_reads {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple platform, coguk_id, run_name, file(bam) from manifest_ch

    output:
    tuple platform, coguk_id, run_name, file(bam), file("${bam}.fasta") into bamfa_manifest_ch

    """
    samtools view ${bam} | awk '{print ">"\$1"\\n"\$10}' > ${bam}.fasta
    """
}

process kraken_bam_reads {
    tag { bam }
    conda "environments/kraken2.yaml"

    input:
    tuple platform, coguk_id, run_name, file(bam), file(bam_fasta) from bamfa_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/k2", pattern: "*k2*", mode: "copy", overwrite: true
    tuple platform, coguk_id, run_name, file(bam), file(bam_fasta), file("${bam_fasta}.k2o.9606.ls") into k2_manifest_ch
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
    tuple platform, coguk_id, run_name, file(bam), file(bam_fasta), file(bam_hum_ls) from k2_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/dh", pattern: "${coguk_id}.${run_name}.dh", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/staging/alignment-clean/", pattern: "${coguk_id}.${run_name}.climb.public.bam", mode: "copy", overwrite: true
    file "${coguk_id}.${run_name}.climb.public.bam"
    file "${coguk_id}.${run_name}.dh" into dh_report_ch

    errorStrategy 'retry' 
    maxRetries 3
    memory { (12 + (2 * task.attempt))+"GB" }

    script:
    if ( platform == "ILLUMINA" )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset sr --bam -o ${coguk_id}.${run_name}.climb.public.bam --trash-minalen 25 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls}
        """
    else if( platform == 'OXFORD_NANOPORE' )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset map-ont --bam -o ${coguk_id}.${run_name}.climb.public.bam --trash-minalen 10 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls}
        """
    else
        error "Invalid alignment mode for technology ${platform}"
}

dh_report_ch
    .collectFile(name: "dehumanised.qc", storeDir: "${params.publish}/staging/dh", keepHeader: true)
