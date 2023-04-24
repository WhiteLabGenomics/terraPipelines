task CollectRNASeqMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File ref_flat
    File ribosomal_intervals


    String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")) + 20
  }

  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      REF_FLAT=~{ref_flat} \
      RIBOSOMAL_INTERVALS= ~{ribosomal_intervals} \
      STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix}.rna_metrics
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File rna_metrics = output_bam_prefix + ".rna_metrics"
  }
}

task CollectMultipleMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index


    String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:3.0.0"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")) + 20
  }

  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    #plots will not be produced if there are no reads
    touch ~{output_bam_prefix}.insert_size_histogram.pdf
    touch ~{output_bam_prefix}.insert_size_metrics
    touch ~{output_bam_prefix}.base_distribution_by_cycle.pdf
    touch ~{output_bam_prefix}.quality_by_cycle.pdf
    touch ~{output_bam_prefix}.quality_distribution.pdf

    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix} \
      REFERENCE_SEQUENCE=~{ref_fasta}
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File alignment_summary_metrics = output_bam_prefix + ".alignment_summary_metrics"
    File insert_size_metrics = output_bam_prefix + ".insert_size_metrics"
    File insert_size_histogram = output_bam_prefix + ".insert_size_histogram.pdf"
    File base_distribution_by_cycle_metrics = output_bam_prefix + ".base_distribution_by_cycle_metrics"
    File base_distribution_by_cycle_pdf = output_bam_prefix + ".base_distribution_by_cycle.pdf"
    File quality_by_cycle_metrics = output_bam_prefix + ".quality_by_cycle_metrics"
    File quality_by_cycle_pdf = output_bam_prefix + ".quality_by_cycle.pdf"
    File quality_distribution_metrics = output_bam_prefix + ".quality_distribution_metrics"
    File quality_distribution_pdf = output_bam_prefix + ".quality_distribution.pdf"
  }
}

