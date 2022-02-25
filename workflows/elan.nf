nextflow.enable.dsl=2

include {save_manifest; resolve_uploads; announce_uploads; samtools_quickcheck; fasta_quickcheck; screen_uploads; rehead_bam; samtools_filter; samtools_index; post_index; samtools_depth; rehead_fasta; swell; post_swell; ocarina_ls} from "../modules/elan.nf"

workflow inbound {
    main:
        if( params.inbound_manifest){
            Channel.fromPath(inbound_manifest)
                .set {inbound_manifest_ch}
        } else {
            save_manifest().out
                .set {inbound_manifest_ch}
        }
        resolve_uploads(inbound_manifest_ch)
        announce_uploads(inbound_manifest_ch, resolve_uploads.out)

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
                "sample_site": it["sourcecode"],             // map sourcecode and sitecode
                "sequencing_site": it["sitecode"]]           //   to less garbage names
            }
            .map { row-> tuple(
                row.coguk_id + ' ' + row.run_name,           // join key
                row,                                         // metadata row
                row.fasta,                                   // input fasta
                row.bam                                      // input bam
	    )}
            .set { manifest_ch }
       
        // Check if the BAM and FASTA are basically garbage
        samtools_quickcheck(manifest_ch)
        fasta_quickcheck(manifest_ch)

        // Screen the uploads for bad fasta and BAM checks taking care to join
        // the outputs of the two independent processes back to the manifest
        // screen_uploads will re-emit the manifest to account for dropped tuples
	screen_uploads(manifest_ch
                                  .join(fasta_quickcheck.out.fasta_rv, by: 0)   // fstatus
                                  .join(samtools_quickcheck.out.bam_rv, by: 0)  // bstatus
        ) | (rehead_bam & rehead_fasta)  // send screened uploads for reheading

        // BAM processing
        rehead_bam.out | samtools_filter | samtools_index | post_index | samtools_depth

        // Join inputs from rehead_fasta, post_index and samtools_depth
        // post_index will drop tuples where indexing failed
        swell(screen_uploads.out.map{ row -> tuple(row[0], row[1]) }     // take the process_key and metadata_row from screen_uploads
                                .join(rehead_fasta.out, by:0)
                                .join(post_index.out.indexed_bam, by:0)  // tuples without an indexed_bam will not be joined
                                .join(samtools_depth.out.bam_depth)
        )

        // Push swell output through post_swell, which will drop tuples that failed swell
        // output can flow right through to ocarina_ls for listing without a join as everything
        // has been nailed down in a tuple through swell
        swell.out.screen_swell | post_swell | ocarina_ls


        // Prepare the Ocarina manifest for elan-ocarina to blast at Majora
        ocarina_ls.out
            .collectFile(name: "ocarina.files.ls", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", sort: false)


        // Collect a QC table that nobody will look at when TQC launches
        swell.out.swell_metrics
            .collectFile(name: "swell.qc.tsv", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", keepHeader: true, sort: false)
        

        // Bung information on all the quickcheck files into a big dumpster for users who wish to go bin diving
        fasta_quickcheck.out.fasta_quickcheck
            .mix( samtools_quickcheck.out.bam_quickcheck )
            .mix( samtools_index.out.bam_index_quickcheck )
            .mix( swell.out.swell_quickcheck )
            .collectFile(name: "elan.quickcheck.ls", storeDir: "${params.artifacts_root}/elan/${params.datestamp}/", sort: false)        
}
