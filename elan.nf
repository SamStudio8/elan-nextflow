#!/usr/bin/env nextflow

params.uploads = "/cephfs/covid/bham/*/upload"
params.publish = "/cephfs/covid/bham/nicholsz/artifacts/elan2"
params.cog_publish = "/cephfs/covid/bham/artifacts/published"
params.minlen = 10000

process save_manifest {
    label 'ocarina'
    conda "environments/ocarina.yaml"

    errorStrategy 'retry'
    maxRetries 1

    output:
    publishDir path: "${params.publish}/staging/summary/${params.datestamp}", pattern: "majora.metadata.tsv", mode: "copy", overwrite: true
    file 'majora.metadata.tsv' into resolve_ch

    """
    ocarina --quiet --env get sequencing --run-name '*' --faster --tsv --task-wait-attempts 75 --task-wait > majora.metadata.tsv
    """
}

process resolve_uploads {

    input:
    file manifest from resolve_ch

    output:
    publishDir path: "${params.publish}/staging/summary/${params.datestamp}", pattern: "files.ls", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/staging/summary/${params.datestamp}", pattern: "files.err", mode: "copy", overwrite: true
    publishDir path: "${params.cog_publish}/elan", pattern: "files.err", mode: "copy", overwrite: true, saveAs: { filename -> "${params.datestamp}.missing.ls" }
    file 'files.ls' into start_ch
    tuple file(manifest), file('files.ls'), file('files.err') into announce_ch
    file 'files.err'
    """
    find ${params.uploads} -type f -name "*fa*" | grep -v '\\.fai\$' | ocarina_resolve.py ${manifest} > files.ls 2> files.err
    """
}

process announce_uploads {
    // https://github.com/COG-UK/dipi-group/issues/89

    label 'ocarina'
    conda "environments/ocarina.yaml"

    validExitStatus 0,1,2,3 // Don't care if this fails, it's just notifying
    //errorStrategy 'retry'
    //maxRetries 1

    input:
    tuple file(manifest), file(q), file(t) from announce_ch

    output:
    publishDir path: "${params.publish}/staging/summary/${params.datestamp}", pattern: "announce.ok", mode: "copy", overwrite: true

    """
    ocarina --env get summary --md > summary.md
    message_uploads.sh ${manifest} ${q} ${t} SHORTSTART \$ELAN_SLACK_HOOK summary.md
    touch announce.ok
    """
}

start_ch
    .splitCsv(header:['is_new', 'coguk_id', 'run_name', 'username', 'pipeuuid', 'autorunname', 'platform', 'dir', 'clabel', 'fasta', 'alabel', 'bam', 'sourcecode', 'sitecode', 'pag', 'tiles', 'adm0', 'adm1_mapped', 'cor_date', 'seq_date'], sep:'\t')
    .filter { row -> row.is_new == "1" }
    .filter { row -> row.fasta.size() > 0 }
    .filter { row -> row.bam.size() > 0 }
    .map { row-> tuple(row.adm0, row.adm1_mapped, row.cor_date, row.seq_date, row.sourcecode, row.sitecode, row.tiles, row.platform, row.pipeuuid, row.username, row.dir, row.run_name, row.coguk_id, file([row.dir, row.fasta].join('/')), file([row.dir, row.bam].join('/'))) }
    .set { manifest_ch }

process samtools_quickcheck {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    validExitStatus 0,1,2,3,4,5,6,7,8

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from manifest_ch

    output:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), env(rv) optional true into validbam_manifest_ch
    file "${coguk_id}.${run_name}.bam.quickcheck" into quickcheck_bam_ch

    shell:
    """
    rv=0
    samtools quickcheck $bam || rv=\$?
    echo "\$rv bam ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.bam.quickcheck
    rv=0
    samtools view $bam > /dev/null || rv=\$?
    echo "\$rv bamv ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" >> ${coguk_id}.${run_name}.bam.quickcheck
    """
}
process fasta_quickcheck {
    tag { fasta }
    label 'bear'

    validExitStatus 0,1,2,3,4,5,6,7,8

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), bstatus from validbam_manifest_ch

    output:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), bstatus, env(rv) optional true into valid_manifest_ch
    file "${coguk_id}.${run_name}.fasta.quickcheck" into quickcheck_fasta_ch

    shell:
    """
    rv=0
    elan_fastacheck.py ${fasta} ${params.minlen} || rv=\$?
    echo "\$rv fasta ${seqsite} ${coguk_id} ${run_name} ${dir}/${fasta}" > ${coguk_id}.${run_name}.fasta.quickcheck
    """
}

process save_uploads {
    tag { fasta }
    label 'bear'

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), bstatus, fstatus from valid_manifest_ch

    output:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file("${coguk_id}.${run_name}.uploaded.fasta"), file("${coguk_id}.${run_name}.uploaded.bam") into validbak_manifest_ch

    errorStrategy 'ignore'

    publishDir path: "${params.publish}/uploaded/alignment", pattern: "${coguk_id}.${run_name}.uploaded.bam", mode: "copy", overwrite: true
    publishDir path: "${params.publish}/uploaded/fasta", pattern: "${coguk_id}.${run_name}.uploaded.fasta", mode: "copy", overwrite: true
    file "${coguk_id}.${run_name}.uploaded.bam"
    file "${coguk_id}.${run_name}.uploaded.fasta"

    script:
    if (fstatus.toInteger() == 0 && bstatus.toInteger() == 0)
        """
        cp ${bam} ${coguk_id}.${run_name}.uploaded.bam
        cp ${fasta} ${coguk_id}.${run_name}.uploaded.fasta
        """
    else
        """
        echo "Cowardly refusing to process ${coguk_id} ${run_name} any further as it has a bad-looking FASTA and/or BAM"
        exit 1
        """

}

process rehead_bam {
    tag { bam }
    label 'bear'
    conda "environments/samtools.yaml"

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from validbak_manifest_ch

    output:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file("${coguk_id}.${run_name}.inbound.bam") into reheaded_manifest_ch

    script:
    """
    samtools view -H --no-PG ${bam} > ${bam}.head
    elan_cleanhead.py ${bam}.head '${coguk_id}.${run_name}' > ${bam}.head.ok
    samtools reheader -P ${bam}.head.ok ${bam} > ${coguk_id}.${run_name}.inbound.bam
    """
}

process samtools_filter {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from reheaded_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/alignment", pattern: "${coguk_id}.${run_name}.climb.bam", mode: "copy", overwrite: true
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file("${coguk_id}.${run_name}.climb.bam") into sorted_manifest_ch

    """
    samtools view -h -F4 ${bam} -o ${coguk_id}.${run_name}.climb.bam
    chmod 644 ${coguk_id}.${run_name}.climb.bam
    """
}

process samtools_index {
    tag { bam }
    label 'bear'
    conda "environments/samtools.yaml"

    errorStrategy 'ignore'

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from sorted_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/alignment", pattern: "${bam.baseName}.bam.bai", mode: "copy", overwrite: true
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), env(rv) into post_index_manifest_ch

    file "${coguk_id}.${run_name}.index.quickcheck" into quickcheck_index_ch
    file "${bam.baseName}.bam.bai" optional true

    script:
    """
    rv=0
    samtools index ${bam} ${bam.baseName}.bam.bai || rv=\$?
    echo "\$rv bam_index ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.index.quickcheck
    """
}
process post_index {
    tag { bam }

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), idx_status from post_index_manifest_ch

    output:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) into indexed_manifest_ch

    errorStrategy 'ignore'

    script:
    if (idx_status.toInteger() == 0)
        """
        echo 'index is good'
        """
    else
        """
        echo "Cowardly refusing to process ${coguk_id} ${run_name} any further as it has a bad-looking BAM"
        exit 1
        """
}

process samtools_depth {
    tag { bam }
    conda "environments/samtools.yaml"
    label 'bear'

    errorStrategy 'retry'
    maxRetries 3
    memory { (3 + (2 * task.attempt))+"GB" }

    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) from indexed_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/depth", pattern: "${coguk_id}.${run_name}.climb.bam.depth", mode: "copy", overwrite: true
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.climb.bam.depth") into swell_manifest_ch

    """
    samtools depth -d0 -a ${bam} > ${bam}.depth
    chmod 644 ${bam}.depth
    """
}

process rehead_fasta {
    label 'bear'
    input:
    tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(depth) from swell_manifest_ch

    output:
    publishDir path : "${params.publish}/staging/fasta/", pattern: "${coguk_id}.${run_name}.climb.fasta", mode: "copy", overwrite: true
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file("${coguk_id}.${run_name}.climb.fasta"), file(bam), file(depth) into swell_ready_manifest_ch

    """
    elan_rehead.py ${fasta} 'COG-UK/${coguk_id}/${seqsite}:${run_name}|${coguk_id}|${adm0}|${adm1}|${sourcesite}|${cor_date}|${seqsite}|${seq_date}' > ${coguk_id}.${run_name}.climb.fasta
    chmod 644 ${coguk_id}.${run_name}.climb.fasta
    """
}

process swell {
    tag { bam }
    conda "environments/swell.yaml"
    label 'bear'

    errorStrategy 'ignore'

    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(depth) from swell_ready_manifest_ch

    output:
    publishDir path: "${params.publish}/staging/qc", pattern: "${coguk_id}.${run_name}.qc", mode: "copy", overwrite: true
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.qc"), env(rv) into postswell_manifest_ch
    file "${coguk_id}.${run_name}.qc" into report_ch
    file "${coguk_id}.${run_name}.swell.quickcheck" into quickcheck_swell_ch

    script:
    if ( tiles == "1" )
        """
        rv=0
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V1/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v1" -x "tileset_reported" "ARTIC-v1" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" --min-pos 1000 --min-pos-allow-total-zero > ${coguk_id}.${run_name}.qc || rv=\$?
        echo "\$rv swell ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.swell.quickcheck
        chmod 644 ${coguk_id}.${run_name}.qc
        """
    else if( tiles == "2" )
        """
        rv=0
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v2" -x "tileset_reported" "ARTIC-v2" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" --min-pos 1000 --min-pos-allow-total-zero > ${coguk_id}.${run_name}.qc || rv=\$?
        echo "\$rv swell ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.swell.quickcheck
        chmod 644 ${coguk_id}.${run_name}.qc
        """
    else if( tiles == "3" )
        """
        rv=0
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v3" -x "tileset_reported" "ARTIC-v3" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" --min-pos 1000 --min-pos-allow-total-zero > ${coguk_id}.${run_name}.qc || rv=\$?
        echo "\$rv swell ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.swell.quickcheck
        chmod 644 ${coguk_id}.${run_name}.qc
        """
    else
        """
        rv=0
        swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --bed "${params.schemegit}/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme.bed" --fasta "${fasta}" -x "tileset_counted" "ARTIC-v2" -x "tileset_reported" "unknown" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" --min-pos 1000 --min-pos-allow-total-zero > ${coguk_id}.${run_name}.qc || rv=\$?
        echo "\$rv swell ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.swell.quickcheck
        chmod 644 ${coguk_id}.${run_name}.qc
        """
}
process post_swell {
    tag { bam }

    input:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(qc), wstatus from postswell_manifest_ch

    output:
    tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(qc) into ocarina_file_manifest_ch

    errorStrategy 'ignore'

    script:
    if (wstatus.toInteger() == 0)
        """
        echo 'swell is good'
        """
    else
        """
        echo "Cowardly refusing to process ${coguk_id} ${run_name} any further as it has a bad-looking BAM"
        exit 1
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

quickcheck_fasta_ch
    .mix( quickcheck_bam_ch )
    .mix( quickcheck_swell_ch )
    .mix( quickcheck_index_ch )
    .collectFile(name: "elan.quickcheck.ls", storeDir: "${params.publish}/staging/summary/${params.datestamp}")

