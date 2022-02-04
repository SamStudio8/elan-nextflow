nextflow.enable.dsl=2

include {save_manifest; resolve_uploads; announce_uploads; samtools_quickcheck; fasta_quickcheck; screen_uploads; rehead_bam; samtools_filter; samtools_index; post_index; samtools_depth; rehead_fasta; swell; post_swell; ocarina_ls} from "../modules/elan.nf"

workflow inbound {
    main:
        save_manifest()
        resolve_uploads(save_manifest.out)
        announce_uploads(save_manifest.out, resolve_uploads.out)

        resolve_uploads.out.elan_files_manifest
            .splitCsv(header:[
                'is_new',
                'coguk_id',
                'run_name',
                'username',
                'pipeuuid',
                'autorunname',
                'platform',
                'dir',
                'clabel',
                'fasta',
                'alabel',
                'bam',
                'sourcecode',
                'sitecode',
                'pag',
                'tiles',
                'adm0',
                'adm1_mapped',
                'cor_date',
                'seq_date'
            ], sep:'\t')
            .filter { row -> row.is_new == "1" }     // only process files marked as new in the manifest
            .filter { row -> row.fasta.size() > 0 }  // with a fasta file name defined
            .filter { row -> row.bam.size() > 0 }    // with a bam file name defined
            .map { it << [
                "fasta": [it["dir"], it["fasta"]].join('/'), // construct the fasta path
                "bam": [it["dir"], it["bam"]].join('/'),     // construct the bam path
                "sample_site": it["sourcesite"],             // map sourcesite and sitecode
                "sequencing_site": it["sitecode"]]           // to less garbage names
            }
            .map { row-> tuple(row, row.fasta, row.bam) }    // create a tuple of metadata_row, fasta and bam
            .set { manifest_ch }
       
        samtools_quickcheck(manifest_ch)
        fasta_quickcheck(manifest_ch)
        screen_uploads(manifest_ch, samtools_quickcheck.out.bam_rv, fasta_quickcheck.out.fasta_rv)
        rehead_bam(manifest_ch, screen_uploads.out.copied_bam)
        samtools_filter(manifest_ch, rehead_bam.out.inbound_bam)
        samtools_index(manifest_ch, samtools_filter.out.filtered_bam)
        post_index(manifest_ch, samtools_index.out.idx_status)
        samtools_depth(manifest_ch, samtools_filter.out.filtered_bam)
        rehead_fasta(manifest_ch, screen_uploads.out.copied_fasta)
        swell(manifest_ch, samtools_filter.out.filtered_bam, rehead_fasta.out.rehead_fasta, samtools_depth.out.bam_depth)
        post_swell(manifest_ch, swell.out.wstatus)
        ocarina_ls(manifest_ch, rehead_fasta.out.rehead_fasta, samtools_filter.out.filtered_bam, swell.out.swell_metrics)

        ocarina_ls.out
            .collectFile(name: "ocarina.files.ls", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", sort: false)

        swell.out.swell_metrics
            .collectFile(name: "swell.qc.tsv", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", keepHeader: true, sort: false)
        
        fasta_quickcheck.out.fasta_quickcheck
            .set { fasta_quickcheck_ch }
        
        samtools_quickcheck.out.bam_quickcheck
            .set { samtools_quickcheck_ch }
        
        swell.out.swell_quickcheck
            .set { swell_quickcheck_ch }
        
        samtools_index.out.bam_index_quickcheck
            .set { bam_index_quickcheck_ch }

        fasta_quickcheck_ch
            .mix( samtools_quickcheck_ch )
            .mix( swell_quickcheck_ch )
            .mix( bam_index_quickcheck_ch )
            .collectFile(name: "elan.quickcheck.ls", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", sort: false)        
}
