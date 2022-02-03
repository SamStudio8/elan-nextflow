process save_manifest {
    label 'ocarina'
    conda "$baseDir/environments/ocarina.yaml"

    errorStrategy 'retry'
    maxRetries 1

    output:
    file 'majora.metadata.tsv'

    publishDir path: "${params.artifacts_root}/elan/${params.datestamp}/", pattern: "majora.metadata.tsv", mode: "copy", overwrite: true

    """
    ocarina --oauth --quiet --env get sequencing --run-name '*' --faster --tsv --task-wait-attempts 75 --task-wait > majora.metadata.tsv
    """
}

process resolve_uploads {

    input:
    file manifest

    output:
    file 'files.ls' 
    file 'files.err'
    
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
    ocarina --oauth --env get summary --md > summary.md
    message_uploads.sh ${manifest} ${files_ls} ${files_err} SHORTSTART \$ELAN_SLACK_HOOK summary.md ${params.datestamp}
    touch announce.ok
    """
}

process samtools_quickcheck {
    tag { bam }
    conda "$baseDir/environments/samtools.yaml"
    label 'bear'

    input:
    val row
    
    output:
    val row, emit: row
    env(rv), emit: bam_rv
    path "${coguk_id}.${run_name}.bam.quickcheck", emit: bam_quickcheck

    shell:
    """
    rv1=0
    samtools quickcheck $bam || rv1=\$?
    echo "\$rv1 bam ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.bam.quickcheck
    rv2=0
    samtools view $bam > /dev/null || rv2=\$?
    echo "\$rv2 bamv ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" >> ${coguk_id}.${run_name}.bam.quickcheck
    rv=\$(( rv1 > rv2 ? rv1 : rv2 ))
    """
}

// process fasta_quickcheck {
//     tag { fasta }
//     label 'bear'

//     input:
//     // tuple adm0, adm1, cor_date, 
//     seq_date, sourcesite
//     seqsite
//     // tiles, platform, pipeuuid, username
//     dir
//     run_name
//     coguk_id
//     file(fasta)
//     // file(bam), bstatus from validbam_manifest_ch

//     output:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), bstatus, 
//     env(rv) optional true
//     file "${coguk_id}.${run_name}.fasta.quickcheck"

//     shell:
//     """
//     rv=0
//     elan_fastacheck.py ${fasta} ${params.minlen} || rv=\$?
//     echo "\$rv fasta ${seqsite} ${coguk_id} ${run_name} ${dir}/${fasta}" > ${coguk_id}.${run_name}.fasta.quickcheck
//     """
// }

// // 2022-01-20 Renamed to screen_uploads as we no longer store the publish/uploaded/ dir as it wastes space
// process screen_uploads {
//     tag { fasta }
//     label 'bear'

//     input:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, 
//     coguk_id
//     file(fasta)
//     file(bam)
//     bstatus
//     fstatus
    
//     output:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id,
//     file("${coguk_id}.${run_name}.uploaded.fasta")
//     file("${coguk_id}.${run_name}.uploaded.bam")

//     errorStrategy 'ignore'

//     // bit pointless now but whatever
//     script:
//     if (fstatus.toInteger() == 0 && bstatus.toInteger() == 0)
//         """
//         cp ${bam} ${coguk_id}.${run_name}.uploaded.bam
//         cp ${fasta} ${coguk_id}.${run_name}.uploaded.fasta
//         """
//     else
//         """
//         echo "Cowardly refusing to process ${coguk_id} ${run_name} any further as it has a bad-looking FASTA and/or BAM"
//         exit 1
//         """

// }

// process rehead_bam {
//     tag { bam }
//     label 'bear'
//     conda "$baseDir/environments/samtools.yaml"

//     input:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir,
//     run_name
//     coguk_id
//     // file(fasta)
//     file(bam)

//     output:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta),
//     file("${coguk_id}.${run_name}.inbound.bam")

//     script:
//     """
//     samtools view -H --no-PG ${bam} > ${bam}.head
//     elan_cleanhead.py ${bam}.head '${coguk_id}.${run_name}' > ${bam}.head.ok
//     samtools reheader -P ${bam}.head.ok ${bam} > ${coguk_id}.${run_name}.inbound.bam
//     """
// }

// process samtools_filter {
//     tag { bam }
//     conda "$baseDir/environments/samtools.yaml"
//     label 'bear'

//     input:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir,
//     run_name
//     coguk_id
//     // file(fasta)
//     file(bam)

//     output:
//     publishDir path: "${params.artifacts_root}/bam/${params.datestamp}/", pattern: "${coguk_id}.${run_name}.climb.bam", mode: "copy", overwrite: true
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta)
//     file("${coguk_id}.${run_name}.climb.bam")

//     """
//     samtools view -h -F4 ${bam} -o ${coguk_id}.${run_name}.climb.bam
//     """
// }

// process samtools_index {
//     tag { bam }
//     label 'bear'
//     conda "$baseDir/environments/samtools.yaml"

//     errorStrategy 'ignore'

//     input:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite,
//     seqsite
//     // tiles, platform, pipeuuid, username,
//     dir
//     run_name
//     coguk_id
//     // file(fasta)
//     file(bam)

//     output:
//     publishDir path: "${params.artifacts_root}/bam/${params.datestamp}/", pattern: "${bam.baseName}.bam.bai", mode: "copy", overwrite: true
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), env(rv) into post_index_manifest_ch

//     file "${coguk_id}.${run_name}.index.quickcheck"
//     file "${bam.baseName}.bam.bai" optional true

//     script:
//     """
//     rv=0
//     samtools index ${bam} ${bam.baseName}.bam.bai || rv=\$?
//     echo "\$rv bam_index ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.index.quickcheck
//     """
// }
// process post_index {
//     tag { bam }

//     input:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, 
//     run_name
//     coguk_id
//     // , file(fasta), file(bam), idx_status from post_index_manifest_ch

//     // output:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam) into indexed_manifest_ch

//     errorStrategy 'ignore'

//     script:
//     if (idx_status.toInteger() == 0)
//         """
//         echo 'index is good'
//         """
//     else
//         """
//         echo "Cowardly refusing to process ${coguk_id} ${run_name} any further as it has a bad-looking BAM"
//         exit 1
//         """
// }

// process samtools_depth {
//     tag { bam }
//     conda "$baseDir/environments/samtools113.yaml"
//     label 'bear'

//     memory '5 GB'

//     input:
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta),
//     file(bam)

//     output:
//     publishDir path: "${params.publish}/staging/depth", pattern: "${coguk_id}.${run_name}.climb.bam.depth", mode: "copy", overwrite: true
//     // tuple adm0, adm1, cor_date, seq_date, sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam),
//     file("${coguk_id}.${run_name}.climb.bam.depth")

//     // 2021-08-25 Updated to use samtools 1.13 for significant improvement to depth algorithm
//     // -d0 removed as depth limit deprecated in samtools 1.13
//     // https://github.com/COG-UK/dipi-group/issues/129
    
//     """
//     samtools depth -a ${bam} > ${bam}.depth
//     """
// }

// process rehead_fasta {
//     label 'bear'
//     input:
//     // tuple adm0, 
//     adm1
//     cor_date
//     seq_date
//     sourcesite
//     seqsite
//     // tiles, platform, pipeuuid, username, dir, 
//     run_name
//     coguk_id
//     file(fasta), file(bam), file(depth) from swell_manifest_ch

//     output:
//     publishDir path : "${params.artifacts_root}/fasta/${params.datestamp}", pattern: "${coguk_id}.${run_name}.climb.fasta", mode: "copy", overwrite: true
//     // tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id,
//     file("${coguk_id}.${run_name}.climb.fasta")
//     // , file(bam), file(depth) into swell_ready_manifest_ch

//     """
//     elan_rehead.py ${fasta} 'COG-UK/${coguk_id}/${seqsite}:${run_name}|${coguk_id}|${adm0}|${adm1}|${sourcesite}|${cor_date}|${seqsite}|${seq_date}' > ${coguk_id}.${run_name}.climb.fasta
//     """
// }


// // Note the allow list for swell uses 'in' rather than exact matching, so NC_045512 will permit NC_045512.2 etc.
// process swell {
//     tag { bam }
//     conda "$baseDir/environments/swell.yaml"
//     label 'bear'

//     errorStrategy 'ignore'

//     input:
//     tuple sourcesite, 
//     seqsite
//     tiles
//     platform
//     // pipeuuid, username
//     dir
//     run_name
//     coguk_id
//     file(fasta)
//     file(bam)
//     file(depth)

//     output:
//     publishDir path: "${params.publish}/staging/qc", pattern: "${coguk_id}.${run_name}.qc", mode: "copy", overwrite: true // still required for ocarina.nf for now, ideally use report_ch version
//     // tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file("${coguk_id}.${run_name}.qc"), env(rv) into postswell_manifest_ch
//     file "${coguk_id}.${run_name}.qc"
//     file "${coguk_id}.${run_name}.swell.quickcheck"

//     // 2022-01-19 Removed dep on artic scheme git, no longer calculating tile depths. TomB will be building tqc -- the next generation QC system. Breaking out smaller subsystems from Majora is the future!
//     // Note that messing with the output here will require a correcting commit to ocarina.nf to ensure the swell output matches the crude readline
//     script:
//     """
//     rv=0
//     swell --ref 'NC_045512' 'NC045512' 'MN908947.3' --depth ${depth} --fasta "${fasta}" -x "tileset_counted" "NA" -x "tileset_reported" "${tiles}" -x "source_site" "${sourcesite}" -x "seq_site" "${seqsite}" -x "platform" "${platform}" -x "datestamp" "${params.datestamp}" --min-pos 1000 --min-pos-allow-total-zero > ${coguk_id}.${run_name}.qc || rv=\$?
//     echo "\$rv swell ${seqsite} ${coguk_id} ${run_name} ${dir}/${bam}" > ${coguk_id}.${run_name}.swell.quickcheck
//     """
// }
// process post_swell {
//     tag { bam }

//     input:
//     // tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir
//     run_name
//     coguk_id
//     // file(fasta), file(bam), file(qc), wstatus from postswell_manifest_ch

//     // output:
//     // tuple sourcesite, seqsite, tiles, platform, pipeuuid, username, dir, run_name, coguk_id, file(fasta), file(bam), file(qc) into ocarina_file_manifest_ch

//     errorStrategy 'ignore'

//     script:
//     if (wstatus.toInteger() == 0)
//         """
//         echo 'swell is good'
//         """
//     else
//         """
//         echo "Cowardly refusing to process ${coguk_id} ${run_name} any further as it has a bad-looking BAM"
//         exit 1
//         """
// }

// // NOTE The entries here need to match the publishDir directives above to make sure Majora knows where the files are
// process ocarina_ls {
//     input:
//     // tuple sourcesite, seqsite, tiles, platform,
//     pipeuuid
//     username
//     // dir,
//     run_name
//     coguk_id
//     file(fasta)
//     file(bam)
//     file(qc)

//     output:
//     file "${coguk_id}.${run_name}.ocarina"
    
//     """
//     echo "${coguk_id}\t${run_name}\t${username}\t${pipeuuid}\t${params.publish}/staging/\t${params.artifacts_root}/fasta/${params.datestamp}/${fasta}\t${params.artifacts_root}/bam/${params.datestamp}/${bam}\tqc/${qc}\t${sourcesite}\t${seqsite}\t${platform}" > ${coguk_id}.${run_name}.ocarina
//     """
// }
