task rebc_sv_consensus_filter_task_1 {
    File input_SV_tsv
    File blacklist_bed
    String id 
    String stub
    String? TALT_thresh1
    String? NALT_thresh1
    String? VAF_thresh1
    String? NALG_thresh1
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    String TALT_thresh="${if defined(TALT_thresh1) then TALT_thresh1 else '4'}"
    String NALT_thresh="${if defined(NALT_thresh1) then NALT_thresh1 else '2'}"
    String VAF_thresh="${if defined(VAF_thresh1) then VAF_thresh1 else '0.1'}"
    String NALG_thresh="${if defined(NALG_thresh1) then NALG_thresh1 else '2'}"
    String ext="${if defined(stub) then '.filtered_SV.tsv' else 'filtered_SV.tsv'}"

    command {
        set -euo pipefail

        echo "${input_SV_tsv}"
        echo "${blacklist_bed}"
        echo "${id}"
        echo "${TALT_thresh}"
        echo "${NALT_thresh}"
        echo "${VAF_thresh}"
        echo "${NALG_thresh}"

        ls -lath
        
        python /opt/src/REBC_SV_filter.py -d ${input_SV_tsv} -b ${blacklist_bed} -i ${id} -s ${stub} -t ${TALT_thresh} -n ${NALT_thresh} -v ${VAF_thresh} -a ${NALG_thresh}

        ls -lath        
    }

    output {
        File consensus_sv_filtered_tsv="${id}.${stub}${ext}"
    }

    runtime {
        docker : "chipstewart/rebc_sv_consensus_filter_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '3'}"
    }

    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }
}

workflow REBC_SV_consensus_filter {

    call rebc_sv_consensus_filter_task_1 

}
