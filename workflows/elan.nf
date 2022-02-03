// ocarina_report_ch
//     .collectFile(name: "ocarina.files.ls", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", sort: false)

// report_ch
//     .collectFile(name: "swell.qc.tsv", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", keepHeader: true, sort: false)

// quickcheck_fasta_ch
//     .mix( quickcheck_bam_ch )
//     .mix( quickcheck_swell_ch )
//     .mix( quickcheck_index_ch )
//     .collectFile(name: "elan.quickcheck.ls", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", sort: false)

nextflow.enable.dsl=2

include {save_manifest; resolve_uploads; announce_uploads; samtools_quickcheck} from "../modules/elan.nf"

workflow inbound {
    main:
        save_manifest()
        resolve_uploads(save_manifest.out)
        announce_uploads(save_manifest.out, resolve_uploads.out)

        resolve_uploads.out.elan_files_manifest
            .splitCsv(header:['is_new', 'coguk_id', 'run_name', 'username', 'pipeuuid', 'autorunname', 'platform', 'dir', 'clabel', 'fasta', 'alabel', 'bam', 'sourcecode', 'sitecode', 'pag', 'tiles', 'adm0', 'adm1_mapped', 'cor_date', 'seq_date'], sep:'\t')
            .filter { row -> row.is_new == "1" }
            .filter { row -> row.fasta.size() > 0 }
            .filter { row -> row.bam.size() > 0 }
            // .map { row-> tuple( file([row.dir, row.fasta].join('/')), file([row.dir, row.bam].join('/'))) }
            .map{ it << ["fasta": file([it["dir"], it["fasta"]].join('/')), "bam": file([it["dir"], it["bam"]].join('/'))] }
            .set { manifest_ch }
        
        samtools_quickcheck(manifest_ch)

}
