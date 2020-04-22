#!/usr/bin/env nextflow

params.uploads = "/cephfs/covid/bham/*/upload"
params.dump = "/cephfs/covid/software/sam/pre-elan/latest.tsv"
params.publish = "/cephfs/covid/bham/nicholsz/artifacts/elan"
params.dhmanifest = "/cephfs/covid/software/sam/dh/20200421/manifest.txt"
params.k2db = "/ramdisk/kraken2db"

process resolve_uploads {

    output:
    file 'files.ls' into start_ch
    """
    find ${params.uploads} -type f -name "*fa*" | grep -v 'fai' | ocarina_resolve.py ${params.dump} > files.ls 2> err
    """
}

start_ch
    .splitCsv(header:['coguk_id', 'run_name', 'username', 'pipeuuid', 'autorunname', 'platform', 'dir', 'clabel', 'fasta', 'alabel', 'bam'], sep:'\t')
    .filter { row -> row.fasta.size() > 0 }
    .filter { row -> row.bam.size() > 0 }
    .map { row-> tuple(row.platform, row.pipeuuid, row.username, row.dir, row.run_name, row.coguk_id, file([row.dir, row.fasta].join('/')), file([row.dir, row.bam].join('/'))) }
    .set { manifest_ch }

process samtools_quickcheck {
    tag { bam }
    conda "environments/samtools.yaml"
    //#label 'bear'

    errorStrategy 'ignore'

    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from manifest_ch

    output:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) into valid_manifest_ch

    """
    samtools quickcheck $bam
    """
}

process save_uploads {
    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from valid_manifest_ch

    output:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) into validbak_manifest_ch

    publishDir path: "${params.publish}/uploaded/alignment", pattern: "${coguk_id}.${run_name}.uploaded.bam", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/uploaded/fasta", pattern: "${coguk_id}.${run_name}.uploaded.fasta", mode: "copy", overwrite: true
    file "${coguk_id}.${run_name}.uploaded.bam"
    file "${coguk_id}.${run_name}.uploaded.fasta"

    """
    cp ${bam} ${coguk_id}.${run_name}.uploaded.bam
    cp ${fasta} ${coguk_id}.${run_name}.uploaded.fasta
    """

}

process extract_bam_reads {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from validbak_manifest_ch

    output:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.bam.fasta") into bamfa_manifest_ch

    """
    samtools view ${bam} | awk '{print ">"\$1"\\n"\$10}' > ${coguk_id}.${run_name}.bam.fasta
    """
}

process kraken_bam_reads {
    tag { bam }
    conda "environments/kraken2.yaml"

    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(bam_fasta) from bamfa_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/k2", pattern: "*k2*", mode: "copy", overwrite: true
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${bam_fasta}.k2o.9606.ls") into k2_manifest_ch
    file "${bam_fasta}.k2o"
    file "${bam_fasta}.k2r"

    cpus 4
    """
    kraken2 --memory-mapping --db ${params.k2db} --threads ${task.cpus} --output ${bam_fasta}.k2o --report ${bam_fasta}.k2r ${bam_fasta} && awk '\$3 == 9606 {print \$2}' ${bam_fasta}.k2o > ${bam_fasta}.k2o.9606.ls
    """
}

process dehumanise_bam {
    tag { bam }
    conda "environments/dehumanizer.yaml"
    label 'bear'

    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(bam_hum_ls) from k2_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/dh", pattern: "${coguk_id}.${run_name}.dh", mode: "copy", overwrite: true
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file("${coguk_id}.${run_name}.dh.bam") into dh_manifest_ch
    file "${coguk_id}.${run_name}.dh" into dh_report_ch

    errorStrategy 'retry' 
    maxRetries 3
    memory { (13 + (2 * task.attempt))+"GB" }

    script:
    if ( platform == "ILL" )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset sr --bam -o ${coguk_id}.${run_name}.dh.bam --trash-minalen 25 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls}
        """
    else if( platform == 'ONT' )
        """
        dehumanise ${params.dhmanifest} ${bam} --preset map-ont --bam -o ${coguk_id}.${run_name}.dh.bam --trash-minalen 10 --log ${coguk_id}.${run_name}.dh --known ${bam_hum_ls}
        """
    else
        error "Invalid alignment mode for technology ${platform}"
}

process samtools_filter_and_sort {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from dh_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/alignment", pattern: "${coguk_id}.${run_name}.climb.bam", mode: "copy", overwrite: true
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file("${coguk_id}.${run_name}.climb.bam") into sorted_manifest_ch

    cpus 4
    memory '5 GB'

    """
    samtools view -h -F4 ${bam} | samtools sort -m1G -@ ${task.cpus} -o ${coguk_id}.${run_name}.climb.bam
    """
}

process samtools_depth {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from sorted_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/depth", pattern: "${coguk_id}.${run_name}.climb.bam.depth", mode: "copy", overwrite: true
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.climb.bam.depth") into swell_manifest_ch

    """
    samtools depth -d0 -a ${bam} > ${bam}.depth
    """
}

process rename_fasta {
    //#label 'bear'
    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(depth) from swell_manifest_ch

    output:
    publishDir path : "${params.publish}/staging/fasta/", pattern: "${coguk_id}.${run_name}.climb.fasta", mode: "copy", overwrite: true
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file("${coguk_id}.${run_name}.climb.fasta"), file(bam), file(depth) into swell_ready_manifest_ch

    """
    cp ${fasta} ${coguk_id}.${run_name}.climb.fasta
    """
}

process swell {
    tag { bam }
    conda "environments/swell.yaml"
    label 'bear'

    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(depth) from swell_ready_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/qc", pattern: "${coguk_id}.${run_name}.qc", mode: "copy", overwrite: true
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.qc") into ocarina_file_manifest_ch
    file "${coguk_id}.${run_name}.qc" into report_ch

    """
    swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme.bed" --fasta "${fasta}" > ${coguk_id}.${run_name}.qc
    """

}

process ocarina_ls {
    input:
    tuple platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(qc) from ocarina_file_manifest_ch

    output:
    file "${coguk_id}.${run_name}.ocarina" into ocarina_report_ch
    
    """
    echo "${coguk_id}\t${run_name}\t${username}\t${pipeuuid}\t${params.publish}/staging/\tfasta/${fasta}\talignment/${bam}\tqc/${qc}" > ${coguk_id}.${run_name}.ocarina
    """
}

ocarina_report_ch
    .collectFile(name: "ocarina.ls", storeDir: "${params.publish}/staging/ocarina")

report_ch
    .collectFile(name: "swell.qc", storeDir: "${params.publish}/staging/swell", keepHeader: true)

dh_report_ch
    .collectFile(name: "dehumanised.qc", storeDir: "${params.publish}/staging/dh", keepHeader: true)
