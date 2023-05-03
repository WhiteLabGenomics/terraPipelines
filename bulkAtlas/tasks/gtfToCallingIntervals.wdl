version 1.0

workflow run_gtfToCallingIntervals {
    call gtfToCallingIntervals
}
######### TASKS #########
task gtfToCallingIntervals {
    input {
        File gtf
        File ref_dict

        String output_name = basename(gtf, ".gtf") + ".exons.interval_list"

        String docker
        String gatk_path
        Int preemptible_count
    }

    String cmd='{print $1 \"\t\" ($2 - 1) \"\t\" $3}'

    command {
        set -e

        echo """args <- commandArgs(trailingOnly = TRUE)
gtf = read.table(args[1], sep='	')
gtf = subset(gtf, V3 == 'exon')
write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), 'exome.bed', quote = F, sep='	', col.names = F, row.names = F)
""" > gtf2exonBed.R

        Rscript gtf2exonBed.R ${gtf}
        
        awk '${cmd}' exome.bed > exome.fixed.bed

        ${gatk_path} \
            BedToIntervalList \
            -I exome.fixed.bed \
            -O ${output_name} \
            -SD ${ref_dict}
    }

    output {
        File interval_list = "${output_name}"
    }

    runtime {
        docker: docker
        preemptible: preemptible_count
    }
}