version 1.0

# https://github.com/labbcb/methseq/blob/master/methseq/workflows/wgbs.wdl

task trimgalore{

    input {
        File fastq_1
        File fastq_2

        Int? five_prime_clip_1
        Int? three_prime_clip_1
        Int? five_prime_clip_2
        Int? three_prime_clip_2

        Int quality = 20

        String cores = "4"
        String? trimgalore_path_override
        String trim_galore_path = select_first([trimgalore_path_override, ""]) + "trim_galore"
    }

    command {
        ~{trim_galore_path} \
            --paired \
            --illumina \
            --gzip \
            ~{'--clip_R1 ' + five_prime_clip_1} \
            ~{'--three_prime_clip_R1 ' + three_prime_clip_1} \
            ~{'--clip_R2 ' + five_prime_clip_2} \
            ~{'--three_prime_clip_R2 ' + three_prime_clip_2} \
            --quality ~{quality} \
            --cores ~{cores} \
            ~{fastq_1} \
            ~{fastq_2}
    }

    output {
        File trim_1 = basename(fastq_1, ".fastq.gz") + "_val_1.fq.gz"
        File trim_2 = basename(fastq_2, ".fastq.gz") + "_val_2.fq.gz"
        File stats_1 = basename(fastq_1) + "_trimming_report.txt"
        File stats_2 = basename(fastq_2) + "_trimming_report.txt"
    }

    runtime {
        docker: "welliton/trimgalore:0.6.5"
        cpu: cores
        memory: "1 GB"
    }
}

workflow trimgalore_workflow {

    input {
        File fastq_1
        File fastq_2

        Int? five_prime_clip_1
        Int? three_prime_clip_1
        Int? five_prime_clip_2
        Int? three_prime_clip_2
        Int quality = 20

        String? trimgalore_path_override
    }

    call trimgalore {
        input:
            fastq_1 = fastq_1,
            fastq_2 = fastq_2,
            five_prime_clip_1 = five_prime_clip_1,
            three_prime_clip_1 = three_prime_clip_1,
            five_prime_clip_2 = five_prime_clip_2,
            three_prime_clip_2 = three_prime_clip_2,
            quality = quality,
            trimgalore_path_override = trimgalore_path_override
    }

    output {
        File trim_1 = trimgalore.trim_1
        File trim_2 = trimgalore.trim_2
        File trim_stats_1 = trimgalore.stats_1
        File trim_stats_2 = trimgalore.stats_2
    }
}