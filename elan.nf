#!/usr/bin/env nextflow

Channel
    .fromPath(params.manifest)
    .splitCsv(header:['dir', 'site', 'coguk_id', 'fasta', 'bam'], sep:'\t')
    .filter { row -> row.fasta.size() > 0 }
    .filter { row -> row.bam.size() > 0 }
    .map { row-> tuple(row.dir, row.site, row.coguk_id, file([row.dir, row.fasta].join('/')), file([row.dir, row.bam].join('/'))) }
    .set { manifest_ch }

process samtools_quickcheck {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    errorStrategy 'ignore'

    input:
    tuple dir, site, coguk_id, file(fasta), file(bam) from manifest_ch

    output:
    tuple dir, site, coguk_id, file(fasta), file(bam) into valid_manifest_ch

    """
    samtools quickcheck $bam
    """
}

process samtools_filter_and_sort {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple dir, site, coguk_id, file(fasta), file(bam) from valid_manifest_ch

    output:
    publishDir path : "${params.publish}/alignment", pattern: "${coguk_id}.climb.bam", mode: "copy"
    tuple dir, site, coguk_id, file(fasta), file("${coguk_id}.climb.bam") into publish_fasta_ch

    cpus 4
    memory '5 GB'

    """
    samtools view -h -F4 ${bam} | samtools sort -m1G -@ ${task.cpus} -o ${coguk_id}.climb.bam
    """
}

process publish_fasta {
    input:
    tuple dir, site, coguk_id, file(fasta), file(bam) from publish_fasta_ch

    output:
    publishDir path : "${params.publish}/fasta", pattern: "${coguk_id}.climb.fasta", mode: "copy"
    tuple dir, site, coguk_id, file("${coguk_id}.climb.fasta"), file(bam) into sorted_manifest_ch

    """
    cp ${fasta} ${coguk_id}.climb.fasta
    """
}

process samtools_depth {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple dir, site, coguk_id, file(fasta), file(bam) from sorted_manifest_ch

    output:
    publishDir path: "${params.publish}/depth", pattern: "${coguk_id}.climb.bam.depth", mode: "copy"
    tuple dir, site, coguk_id, file(fasta), file(bam), file("${coguk_id}.climb.bam.depth") into swell_manifest_ch

    """
    samtools depth -d0 -a ${bam} > ${bam}.depth
    """
}

process swell {
    tag { bam }
    conda "environments/swell.yaml"
    label 'bear'

    input:
    tuple dir, site, coguk_id, file(fasta), file(bam), file(depth) from swell_manifest_ch

    output:
    publishDir path: "${params.publish}/qc", pattern: "${coguk_id}.qc", mode: "copy"
    file "${coguk_id}.qc" into report_ch

    """
    swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme.bed" --fasta "${fasta}" > ${coguk_id}.qc
    """

}

report_ch
    .collectFile(name: "test.qc", storeDir: "${params.publish}/qc", keepHeader: true)
