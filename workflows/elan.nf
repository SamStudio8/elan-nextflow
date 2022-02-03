// start_ch
//     .splitCsv(header:['is_new', 'coguk_id', 'run_name', 'username', 'pipeuuid', 'autorunname', 'platform', 'dir', 'clabel', 'fasta', 'alabel', 'bam', 'sourcecode', 'sitecode', 'pag', 'tiles', 'adm0', 'adm1_mapped', 'cor_date', 'seq_date'], sep:'\t')
//     .filter { row -> row.is_new == "1" }
//     .filter { row -> row.fasta.size() > 0 }
//     .filter { row -> row.bam.size() > 0 }
//     .map { row-> tuple(row.adm0, row.adm1_mapped, row.cor_date, row.seq_date, row.sourcecode, row.sitecode, row.tiles, row.platform, row.pipeuuid, row.username, row.dir, row.run_name, row.coguk_id, file([row.dir, row.fasta].join('/')), file([row.dir, row.bam].join('/'))) }
//     .set { manifest_ch }


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

include {save_manifest; resolve_uploads} from "../modules/elan.nf"

workflow inbound {
    main:
        save_manifest()
        resolve_uploads(save_manifest.out)
}
