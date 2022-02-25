process save_manifest {
    label 'ocarina'
    conda "$baseDir/environments/ocarina.yaml"

    errorStrategy 'retry'
    maxRetries 1

    output:
    file 'majora.metadata.tsv'

    publishDir path: "${params.artifacts_root}/elan/${params.datestamp}/", pattern: "majora.metadata.tsv", mode: "copy", overwrite: true

    """
    ocarina --oauth --quiet --profile ${params.ocarina_profile} get sequencing --run-name '*' --faster --tsv --task-wait-attempts 90 --task-wait > majora.metadata.tsv
    """
}

process resolve_uploads {

    input:
    file manifest

    output:
    path 'files.ls', emit: elan_files_manifest
    path 'files.err', emit: elan_files_unmatched
    
    publishDir path: "${params.artifacts_root}/elan/${params.datestamp}/", pattern: "files.ls", mode: "copy", overwrite: true, saveAs: { filename -> "elan.manifest.ls" }
    publishDir path: "${params.artifacts_root}/elan/${params.datestamp}/", pattern: "files.err", mode: "copy", overwrite: true, saveAs: { filename -> "elan.missing.ls" }
        
    
    """
    find ${params.uploads} -type f -name "*fa*" | grep -v '\\.fai\$' | ocarina_resolve.py --metadata ${manifest} --user-field ${params.uploads_usern} > files.ls 2> files.err
    """
}

process announce_uploads {
    // https://github.com/COG-UK/dipi-group/issues/89

    label 'ocarina'
    conda "$baseDir/environments/ocarina.yaml"

    errorStrategy 'ignore'

    input:
    file(manifest)
    file(files_ls)
    file(files_err)

    publishDir path: "${params.artifacts_root}/elan/${params.datestamp}", pattern: "announce.ok", mode: "copy", overwrite: true

    """
    ocarina --oauth --profile ${params.ocarina_profile} get summary --md > summary.md
    message_uploads.sh ${manifest} ${files_ls} ${files_err} SHORTSTART \$ELAN_SLACK_HOOK summary.md ${params.datestamp}
    touch announce.ok
    """
}

process samtools_quickcheck {
    tag { process_key }
    conda "$baseDir/environments/samtools.yaml"
    label 'bear'

    input:
    tuple val(process_key), val(row), path(fasta), path(bam)
    
    output:
    tuple val(process_key), env(rv), emit: bam_rv
    path "${row.coguk_id}.${row.run_name}.bam.quickcheck", emit: bam_quickcheck

    shell:
    """
    rv1=0
    samtools quickcheck ${bam} || rv1=\$?
    echo "\$rv1 bam ${row.sequencing_site} ${row.coguk_id} ${row.run_name} ${row.dir}/${bam}" > ${row.coguk_id}.${row.run_name}.bam.quickcheck
    rv2=0
    samtools view ${bam} > /dev/null || rv2=\$?
    echo "\$rv2 bamv ${row.sequencing_site} ${row.coguk_id} ${row.run_name} ${row.dir}/${bam}" >> ${row.coguk_id}.${row.run_name}.bam.quickcheck
    rv=\$(( rv1 > rv2 ? rv1 : rv2 ))
    """
}

// NOTE must cover any potential exit conditions of rehead_fasta as that process is not allowed to exit non-zero
process fasta_quickcheck {
    tag { process_key }
    label 'bear'

    input:
    tuple val(process_key), val(row), path(fasta), path(bam)

    output:
    tuple val(process_key), env(rv), emit: fasta_rv
    path "${row.coguk_id}.${row.run_name}.fasta.quickcheck", emit: fasta_quickcheck

    shell:
    """
    rv=0
    elan_fastacheck.py ${fasta} ${params.minlen} || rv=\$?
    echo "\$rv fasta ${row.sequencing_site} ${row.coguk_id} ${row.run_name} ${row.dir}/${fasta}" > ${row.coguk_id}.${row.run_name}.fasta.quickcheck
    """
}

// 2022-01-20 Renamed to screen_uploads as we no longer store the publish/uploaded/ dir as it wastes space
process screen_uploads {
    tag { process_key }
    label 'bear'

    input:
    tuple val(process_key), val(row), path(fasta), path(bam), val(fstatus), val(bstatus)
    
    output:
    tuple val(process_key), val(row), path("${row.coguk_id}.${row.run_name}.uploaded.fasta"), path("${row.coguk_id}.${row.run_name}.uploaded.bam"), emit: screened_uploads

    errorStrategy 'ignore'

    // copy is a bit pointless now we dont store raw backup but whatever
    script:
    if (fstatus.toInteger() == 0 && bstatus.toInteger() == 0)
        """
        cp ${bam} ${row.coguk_id}.${row.run_name}.uploaded.bam
        cp ${fasta} ${row.coguk_id}.${row.run_name}.uploaded.fasta
        """
    else
        """
        echo "Cowardly refusing to process ${row.coguk_id} ${row.run_name} any further as it has a bad-looking FASTA and/or BAM"
        exit 1
        """

}

process rehead_bam {
    tag { process_key }
    label 'bear'
    conda "$baseDir/environments/samtools.yaml"

    input:
    tuple val(process_key), val(row), path(copied_fasta), path(copied_bam)

    output:
    tuple val(process_key), val(row), path("${row.coguk_id}.${row.run_name}.inbound.bam"), emit: inbound_bam

    script:
    """
    samtools view -H --no-PG ${copied_bam} > ${copied_bam}.head
    elan_cleanhead.py ${copied_bam}.head '${row.coguk_id}.${row.run_name}' > ${copied_bam}.head.ok
    samtools reheader -P ${copied_bam}.head.ok ${copied_bam} > ${row.coguk_id}.${row.run_name}.inbound.bam
    """
}

process samtools_filter {
    tag { process_key }
    conda "$baseDir/environments/samtools.yaml"
    label 'bear'

    input:
    tuple val(process_key), val(row), path(inbound_bam)

    output:
    tuple val(process_key), val(row), path("${row.coguk_id}.${row.run_name}.climb.bam"), emit: filtered_bam

    publishDir path: "${params.artifacts_root}/bam/${params.datestamp}/", pattern: "${row.coguk_id}.${row.run_name}.climb.bam", mode: "copy", overwrite: true

    """
    samtools view -h -F4 ${inbound_bam} -o ${row.coguk_id}.${row.run_name}.climb.bam
    """
}

process samtools_index {
    tag { process_key }
    label 'bear'
    conda "$baseDir/environments/samtools.yaml"

    errorStrategy 'ignore'

    input:
    tuple val(process_key), val(row), path(filtered_bam)

    output:
    val(row), emit: metadata
    tuple val(process_key), path(filtered_bam), emit: indexed_bam
    env(rv), emit: idx_status
    path "${row.coguk_id}.${row.run_name}.index.quickcheck", emit: bam_index_quickcheck
    path "${filtered_bam.baseName}.bam.bai", emit: bam_bai optional true

    publishDir path: "${params.artifacts_root}/bam/${params.datestamp}/", pattern: "${filtered_bam.baseName}.bam.bai", mode: "copy", overwrite: true

    // todo sw has changed these but the quickcheck paths need to point to user bam
    script: //bam.bai specified so that it is created within working dir rather than another, separate working dir
    """
    rv=0
    samtools index ${filtered_bam} ${filtered_bam.baseName}.bam.bai || rv=\$?
    echo "\$rv bam_index ${row.sequencing_site} ${row.coguk_id} ${row.run_name} ${row.dir}/${filtered_bam}" > ${row.coguk_id}.${row.run_name}.index.quickcheck
    """
}

process post_index {
    tag { process_key }

    input:
    val(row)
    tuple val(process_key), path(indexed_bam)
    val(idx_status)
    path(quickcheck)
    path(index)

    output:
    val(row), emit: metadata
    tuple val(process_key), path(indexed_bam), emit: indexed_bam

    errorStrategy 'ignore'

    script:
    if (idx_status.toInteger() == 0)
        """
        echo 'index is good'
        """
    else
        """
        echo "Cowardly refusing to process ${row.coguk_id} ${row.run_name} any further as it has a bad-looking BAM"
        exit 1
        """
}

process samtools_depth {
    tag { process_key }
    conda "$baseDir/environments/samtools113.yaml"
    label 'bear'

    memory '5 GB'

    input:
    val(row)
    tuple val(process_key), path(indexed_bam)

    output:
    tuple val(process_key), path("${row.coguk_id}.${row.run_name}.climb.bam.depth"), emit: bam_depth

    // 2021-08-25 Updated to use samtools 1.13 for significant improvement to depth algorithm
    // -d0 removed as depth limit deprecated in samtools 1.13
    // https://github.com/COG-UK/dipi-group/issues/129
    
    """
    samtools depth -a ${indexed_bam} > ${indexed_bam}.depth
    """
}

process rehead_fasta {
    tag { process_key }
    label 'bear'

    input:
    tuple val(process_key), val(row), path(copied_fasta), path(copied_bam)
    
    output:
    tuple val(process_key), path("${row.coguk_id}.${row.run_name}.climb.fasta"), emit: rehead_fasta

    publishDir path : "${params.artifacts_root}/fasta/${params.datestamp}", pattern: "${row.coguk_id}.${row.run_name}.climb.fasta", mode: "copy", overwrite: true

    """
    elan_rehead.py ${copied_fasta} 'COG-UK/${row.coguk_id}/${row.sequencing_site}:${row.run_name}|${row.coguk_id}|${row.adm0}|${row.adm1_mapped}|${row.sample_site}|${row.cor_date}|${row.sequencing_site}|${row.seq_date}' > ${row.coguk_id}.${row.run_name}.climb.fasta
    """
}


// Note the allow list for swell uses 'in' rather than exact matching, so NC_045512 will permit NC_045512.2 etc.
process swell {
    tag { process_key }
    conda "$baseDir/environments/swell.yaml"
    label 'bear'

    errorStrategy 'ignore'

    input:
    tuple val(process_key), val(row), path(rehead_fasta), path(filtered_bam), path(depth)

    output:
    path "${row.coguk_id}.${row.run_name}.qc", emit: swell_metrics
    path "${row.coguk_id}.${row.run_name}.swell.quickcheck", emit: swell_quickcheck
    tuple val(process_key), val(row), path(rehead_fasta), path(filtered_bam), path("${row.coguk_id}.${row.run_name}.qc"), env(rv), emit: screen_swell

    publishDir path: "${params.publish}/staging/qc", pattern: "${row.coguk_id}.${row.run_name}.qc", mode: "copy", overwrite: true // still required for ocarina.nf for now, ideally use report_ch version

    // 2022-01-19 Removed dep on artic scheme git, no longer calculating tile depths. TomB will be building tqc -- the next generation QC system. Breaking out smaller subsystems from Majora is the future!
    // Note that messing with the output here will require a correcting commit to ocarina.nf to ensure the swell output matches the crude readline
    script:
    """
    rv=0
    swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --fasta "${rehead_fasta}" -x "tileset_counted" "NA" -x "tileset_reported" "NA" -x "source_site" "${row.sample_site}" -x "seq_site" "${row.sequencing_site}" -x "platform" "${row.platform}" -x "datestamp" "${params.datestamp}" --min-pos 1000 --min-pos-allow-total-zero > ${row.coguk_id}.${row.run_name}.qc || rv=\$?
    echo "\$rv swell ${row.sequencing_site} ${row.coguk_id} ${row.run_name} ${row.dir}/${filtered_bam}" > ${row.coguk_id}.${row.run_name}.swell.quickcheck
    """
}

process post_swell {
    tag { process_key }

    input:
    tuple val(process_key), val(row), path(reheaded_fasta), path(indexed_bam), path(swell_metrics), val(wstatus)

    errorStrategy 'ignore'

    output:
    tuple val(process_key), val(row), path(reheaded_fasta), path(indexed_bam), path(swell_metrics)

    script:
    if (wstatus.toInteger() == 0)
        """
        echo 'swell is good'
        """
    else
        """
        echo "Cowardly refusing to process ${row.coguk_id} ${row.run_name} any further as it has a bad-looking BAM"
        exit 1
        """
}

// NOTE The entries here need to match the publishDir directives above to make sure Majora knows where the files are
process ocarina_ls {
    tag { process_key }
    
    input:
    tuple val(process_key), val(row), path(rehead_fasta), path(filtered_bam), path(swell_metrics)

    output:
    path "${row.coguk_id}.${row.run_name}.ocarina"
    
    """
    echo "${row.coguk_id}\t${row.run_name}\t${row.username}\t${row.pipeuuid}\t${params.publish}/staging/\t${params.artifacts_root}/fasta/${params.datestamp}/${rehead_fasta}\t${params.artifacts_root}/bam/${params.datestamp}/${filtered_bam}\tqc/${swell_metrics}\t${row.sample_site}\t${row.sequencing_site}\t${row.platform}" > ${row.coguk_id}.${row.run_name}.ocarina
    """
}
