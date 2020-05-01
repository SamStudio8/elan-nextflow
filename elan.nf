#!/usr/bin/env nextflow

params.uploads = "/cephfs/covid/bham/*/upload"
params.dump = "/cephfs/covid/software/sam/pre-elan/latest.tsv"
params.publish = "/cephfs/covid/bham/nicholsz/artifacts/elan2"

process resolve_uploads {

    output:
    publishDir path: "${params.publish}/staging/summary/${params.datestamp}", pattern: "files.ls", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/staging/summary/${params.datestamp}", pattern: "files.err", mode: "copy", overwrite: true
    file 'files.ls' into start_ch
    file 'files.err'
    """
    find ${params.uploads} -type f -name "*fa*" | grep -v 'fai' | ocarina_resolve.py ${params.dump} > files.ls 2> files.err
    """
}

start_ch
    .splitCsv(header:['is_new', 'coguk_id', 'run_name', 'username', 'pipeuuid', 'autorunname', 'platform', 'dir', 'clabel', 'fasta', 'alabel', 'bam', 'sourcecode', 'sitecode', 'pag', 'tiles'], sep:'\t')
    .filter { row -> row.is_new == "1" }
    .filter { row -> row.fasta.size() > 0 }
    .filter { row -> row.bam.size() > 0 }
    .map { row-> tuple(row.sourcecode, row.sitecode, row.tiles, row.platform, row.pipeuuid, row.username, row.dir, row.run_name, row.coguk_id, file([row.dir, row.fasta].join('/')), file([row.dir, row.bam].join('/'))) }
    .set { manifest_ch }

process samtools_quickcheck {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    errorStrategy 'ignore'

    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from manifest_ch

    output:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) into valid_manifest_ch

    """
    samtools quickcheck $bam
    """
}

process save_uploads {
    label 'bear'

    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from valid_manifest_ch

    output:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) into validbak_manifest_ch

    publishDir path: "${params.publish}/uploaded/alignment", pattern: "${coguk_id}.${run_name}.uploaded.bam", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/uploaded/fasta", pattern: "${coguk_id}.${run_name}.uploaded.fasta", mode: "copy", overwrite: true
    file "${coguk_id}.${run_name}.uploaded.bam"
    file "${coguk_id}.${run_name}.uploaded.fasta"

    """
    cp ${bam} ${coguk_id}.${run_name}.uploaded.bam
    cp ${fasta} ${coguk_id}.${run_name}.uploaded.fasta
    """

}

process samtools_filter_and_sort {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from validbak_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/alignment", pattern: "${coguk_id}.${run_name}.climb.bam", mode: "copy", overwrite: true
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file("${coguk_id}.${run_name}.climb.bam") into sorted_manifest_ch

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
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from sorted_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/depth", pattern: "${coguk_id}.${run_name}.climb.bam.depth", mode: "copy", overwrite: true
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.climb.bam.depth") into swell_manifest_ch

    """
    samtools depth -d0 -a ${bam} > ${bam}.depth
    """
}

process rename_fasta {
    label 'bear'
    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(depth) from swell_manifest_ch

    output:
    publishDir path : "${params.publish}/staging/fasta/", pattern: "${coguk_id}.${run_name}.climb.fasta", mode: "copy", overwrite: true
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file("${coguk_id}.${run_name}.climb.fasta"), file(bam), file(depth) into swell_ready_manifest_ch

    """
    cp ${fasta} ${coguk_id}.${run_name}.climb.fasta
    """
}

process swell {
    tag { bam }
    conda "environments/swell.yaml"
    label 'bear'

    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(depth) from swell_ready_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/qc", pattern: "${coguk_id}.${run_name}.qc", mode: "copy", overwrite: true
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.qc") into ocarina_file_manifest_ch
    file "${coguk_id}.${run_name}.qc" into report_ch

    script:
    if ( tiles == "1" )
        """
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V1/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v1" -x "tileset_reported" "ARTIC-v1" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" > ${coguk_id}.${run_name}.qc
        """
    else if( tiles == "2" )
        """
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v2" -x "tileset_reported" "ARTIC-v2" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" > ${coguk_id}.${run_name}.qc
        """
    else if( tiles == "3" )
        """
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v3" -x "tileset_reported" "ARTIC-v3" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" > ${coguk_id}.${run_name}.qc
        """
    else
        """
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v2" -x "tileset_reported" "unknown" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" > ${coguk_id}.${run_name}.qc
        """
}

process ocarina_ls {
    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(qc) from ocarina_file_manifest_ch

    output:
    file "${coguk_id}.${run_name}.ocarina" into ocarina_report_ch
    
    """
    echo "${coguk_id}\t${run_name}\t${username}\t${pipeuuid}\t${params.publish}/staging/\tfasta/${fasta}\talignment/${bam}\tqc/${qc}\t${sourcesite}\t${seqsite}\t${platform}" > ${coguk_id}.${run_name}.ocarina
    """
}

ocarina_report_ch
    .collectFile(name: "ocarina.files.ls", storeDir: "${params.publish}/staging/summary/${params.datestamp}")

report_ch
    .collectFile(name: "swell.qc.tsv", storeDir: "${params.publish}/staging/summary/${params.datestamp}", keepHeader: true)

