#!/usr/bin/env nextflow

// params.manifest = "/cephfs/covid/bham/nicholsz/artifacts/elan2/staging/summary/null/ocarina.files.ls"

Channel
    .fromPath(params.manifest)
    .splitCsv(header:['coguk_id', 'run_name', 'username', 'pipeuuid', 'dir', 'fasta', 'bam', 'qc', 'sourcesite', 'seqsite', 'platform'], sep:'\t')
    .map { row-> tuple(row.coguk_id, row.run_name, row.username, row.pipeuuid, [row.dir, row.fasta].join('/'), [row.dir, row.bam].join('/'), [row.dir, row.qc].join('/'), row.sourcesite, row.seqsite, row.platform) }
    .set { manifest_ch }

process play_ocarina {
    tag { bam }
    label 'ocarina'
    conda "environments/ocarina.yaml"

    input:
    tuple coguk_id, run_name, username, pipeuuid, fasta, bam, qc, sourcesite, seqsite, platform from manifest_ch

    cpus 8 //# massively over-request local cores to prevent sending too much to API at once

    script:
    """
    ocarina --angry --env put file --path "${fasta}" --type consensus --i-have-bad-files --source-artifact "sequencing-dummy-reads-${run_name}" --bridge-artifact "${coguk_id}" --pipeline "bioinfo-${run_name}" Pipeline unknown 0 --sudo-as "${username}" --full-path --node climb --publish-group "COG-UK/${coguk_id}/${seqsite}:${run_name}" --meta sequencing platform "${platform}";
    ocarina --angry --env put file --path "${bam}" --type alignment --i-have-bad-files --source-artifact "sequencing-dummy-reads-${run_name}" --bridge-artifact "${coguk_id}" --pipeline "bioinfo-${run_name}" Pipeline unknown 0 --sudo-as "${username}" --full-path --node climb --publish-group "COG-UK/${coguk_id}/${seqsite}:${run_name}";

	tail -n1 ${qc} | while read fp num_seqs num_bases pc_acgt pc_masked pc_invalid longest_gap longest_ungap bp num_pos mean_cov pc_pos_cov_gte1 pc_pos_cov_gte5 pc_pos_cov_gte10 pc_pos_cov_gte20 pc_pos_cov_gte50 pc_pos_cov_gte100 pc_pos_cov_gte200 pc_tiles_gte1 pc_tiles_gte5 pc_tiles_gte10 pc_tiles_gte20 pc_tiles_gte50 pc_tiles_gte100 pc_tiles_gte200 n_tiles tile_vector tileset_counted tileset_reported; \n
	do
	   ocarina --env put metric --artifact-path ${fasta} \
                                                     -m sequence num_seqs \$num_seqs \
                                                     -m sequence num_bases \$num_bases \
                                                     -m sequence pc_acgt \$pc_acgt \
                                                     -m sequence pc_masked \$pc_masked \
                                                     -m sequence pc_invalid \$pc_invalid \
                                                     -m sequence longest_gap \$longest_gap \
                                                     -m sequence longest_ungap \$longest_ungap;
	   ocarina --env put metric --artifact-path ${bam} \
                                                     -m mapping num_pos \$num_pos \
                                                     -m mapping mean_cov \$mean_cov \
                                                     -m mapping pc_pos_cov_gte1 \$pc_pos_cov_gte1 \
                                                     -m mapping pc_pos_cov_gte5 \$pc_pos_cov_gte5 \
                                                     -m mapping pc_pos_cov_gte10 \$pc_pos_cov_gte10 \
                                                     -m mapping pc_pos_cov_gte20 \$pc_pos_cov_gte20 \
                                                     -m mapping pc_pos_cov_gte50 \$pc_pos_cov_gte50 \
                                                     -m mapping pc_pos_cov_gte100 \$pc_pos_cov_gte100 \
                                                     -m mapping pc_pos_cov_gte200 \$pc_pos_cov_gte200;
                                                     #-m tile-mapping n_tiles \$n_tiles \
                                                     #-m tile-mapping tile_set \$tileset_reported \
                                                     #-m tile-mapping pc_tiles_medcov_gte1 \$pc_tiles_gte1 \
                                                     #-m tile-mapping pc_tiles_medcov_gte5 \$pc_tiles_gte5 \
                                                     #-m tile-mapping pc_tiles_medcov_gte10 \$pc_tiles_gte10 \
                                                     #-m tile-mapping pc_tiles_medcov_gte20 \$pc_tiles_gte20 \
                                                     #-m tile-mapping pc_tiles_medcov_gte50 \$pc_tiles_gte50 \
                                                     #-m tile-mapping pc_tiles_medcov_gte100 \$pc_tiles_gte100 \
                                                     #-m tile-mapping pc_tiles_medcov_gte200 \$pc_tiles_gte200;
        done
        ocarina --angry --env put qc --publish-group "COG-UK/${coguk_id}/${seqsite}:${run_name}" --test-name 'cog-uk-elan-minimal-qc' --test-version 1;
        ocarina --angry --env put qc --publish-group "COG-UK/${coguk_id}/${seqsite}:${run_name}" --test-name 'cog-uk-high-quality-public' --test-version 1;
    """
}
