process play_ocarina {
    tag { bam }
    label 'ocarina'
    conda "$baseDir/environments/ocarina.yaml"

    maxRetries 2
    errorStrategy { sleep(Math.pow(2, task.attempt) * 500 as long); task.attempt <= maxRetries ? 'retry' : 'terminate' }

    input:
    val row

    // note that --node climb requires majora to have a corresponding DigitalResourceNode with unique_name set to climb
    // additionally cog-uk-elan-minimal-qc and cog-uk-high-quality-public must be loaded with python manage.py load_qctest
    script:
    """
    ocarina --oauth --angry --env put file --path "${row.fasta}" --type consensus --i-have-bad-files --source-artifact "sequencing-dummy-reads-${row.run_name}" --bridge-artifact "${row.coguk_id}" --pipeline-hook "bioinfo-${row.run_name}" --sudo-as "${row.username}" --full-path --node climb --publish-group "COG-UK/${row.coguk_id}/${row.seqsite}:${row.run_name}" --meta sequencing platform "${row.platform}";
    ocarina --oauth --angry --env put file --path "${row.bam}" --type alignment --i-have-bad-files --source-artifact "sequencing-dummy-reads-${row.run_name}" --bridge-artifact "${row.coguk_id}" --pipeline-hook "bioinfo-${row.run_name}" --sudo-as "${row.username}" --full-path --node climb --publish-group "COG-UK/${row.coguk_id}/${row.seqsite}:${row.run_name}";

	tail -n1 ${row.qc} | while read fasfp num_seqs num_bases pc_acgt pc_masked pc_invalid longest_gap longest_ungap bamfp num_pos mean_cov pc_pos_cov_gte1 pc_pos_cov_gte5 pc_pos_cov_gte10 pc_pos_cov_gte20 pc_pos_cov_gte50 pc_pos_cov_gte100 pc_pos_cov_gte200 pc_tiles_gte1 pc_tiles_gte5 pc_tiles_gte10 pc_tiles_gte20 pc_tiles_gte50 pc_tiles_gte100 pc_tiles_gte200 n_tiles tile_vector tileset_counted tileset_reported; \n
	do
	   ocarina --oauth --env put metric --artifact-path ${row.fasta} \
                                                     -m sequence num_seqs \$num_seqs \
                                                     -m sequence num_bases \$num_bases \
                                                     -m sequence pc_acgt \$pc_acgt \
                                                     -m sequence pc_masked \$pc_masked \
                                                     -m sequence pc_invalid \$pc_invalid \
                                                     -m sequence longest_gap \$longest_gap \
                                                     -m sequence longest_ungap \$longest_ungap;
	   ocarina --oauth --env put metric --artifact-path ${row.bam} \
                                                     -m mapping num_pos \$num_pos \
                                                     -m mapping mean_cov \$mean_cov \
                                                     -m mapping pc_pos_cov_gte1 \$pc_pos_cov_gte1 \
                                                     -m mapping pc_pos_cov_gte5 \$pc_pos_cov_gte5 \
                                                     -m mapping pc_pos_cov_gte10 \$pc_pos_cov_gte10 \
                                                     -m mapping pc_pos_cov_gte20 \$pc_pos_cov_gte20 \
                                                     -m mapping pc_pos_cov_gte50 \$pc_pos_cov_gte50 \
                                                     -m mapping pc_pos_cov_gte100 \$pc_pos_cov_gte100 \
                                                     -m mapping pc_pos_cov_gte200 \$pc_pos_cov_gte200;
        done
        ocarina --oauth --angry --env put qc --publish-group "COG-UK/${row.coguk_id}/${row.seqsite}:${row.run_name}" --test-name 'cog-uk-elan-minimal-qc' --test-version 1;
        ocarina --oauth --angry --env put qc --publish-group "COG-UK/${row.coguk_id}/${row.seqsite}:${row.run_name}" --test-name 'cog-uk-high-quality-public' --test-version 1;
    """
}
