version 1.0

import "fastqc.wdl" as fastqc_v1
import "trimgalore.wdl" as trimgalore_v1
import "star.wdl" as star_v1
import "rnaseqc2.wdl" as rnaseqc2_v1
import "rsem.wdl" as rsem_v1


workflow RNA_preprocessing_pipeline {

  input {

    File fastq1
    File fastq2
    String sample_id
    String? trimgalore_path_override ##TODO

    #rannotation GTF
    File genes_gtf="gs://ccle_default_params/references_gtex_gencode.v29.GRCh38.ERCC.genes.collapsed_only.gtf"

    #star_v1 index
    File star_index="gs://ccle_default_params/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v29_oh100.tar.gz"

    #rsem index
    File rsem_reference="gs://ccle_default_params/rsem_reference_GRCh38_gencode29_ercc.tar.gz"

  }


  call fastqc_v1.fastqc as raw_fastqc1 {
    input:
      seqFile=fastq1,
      outdirPath="."
  }

  call fastqc_v1.fastqc as raw_fastqc2 {
    input:
      seqFile=fastq2,
      outdirPath="."
  }

  call trimgalore_v1.trimgalore as trimgalore {
      input:
          fastq_1 = fastq1,
          fastq_2 = fastq2,
          trimgalore_path_override = trimgalore_path_override ## TODO
    }

  call fastqc_v1.fastqc as cleaned_fastqc1 {
    input:
      seqFile=trimgalore.trim_1,
      outdirPath="."
  }

  call fastqc_v1.fastqc as cleaned_fastqc2 {
    input:
      seqFile=trimgalore.trim_2,
      outdirPath="."
  }

  call star_v1.star as star {
    input:
      prefix=sample_id,
      fastq1=trimgalore.trim_1,
      fastq2=trimgalore.trim_2,
      star_index=star_index
  }

  call rnaseqc2_v1.rnaseqc2 as rnaseqc2 {
    input:
      bam_file=star.bam_file,
      genes_gtf=genes_gtf,
      sample_id=sample_id
  }

  call rsem_v1.rsem as rsem {
    input:
      transcriptome_bam=star.transcriptome_bam,
      prefix=sample_id,
      rsem_reference=rsem_reference,
  }


  output {
    #fastqc raw data
    File raw_htmlReport1 = raw_fastqc1.htmlReport
    File raw_reportZip1 = raw_fastqc1.reportZip
    File raw_htmlReport2 = raw_fastqc2.htmlReport
    File raw_reportZip2 = raw_fastqc2.reportZip

    #trimgalore
    File trim_1 = trimgalore.trim_1
    File trim_2 = trimgalore.trim_2
    File trim_stats_1 = trimgalore.stats_1
    File trim_stats_2 = trimgalore.stats_2

    #fastqc cleaned data
    File cleaned_htmlReport1 = cleaned_fastqc1.htmlReport
    File cleaned_reportZip1 = cleaned_fastqc1.reportZip
    File cleaned_htmlReport2 = cleaned_fastqc2.htmlReport
    File cleaned_reportZip2 = cleaned_fastqc2.reportZip

    #star
    File bam_file=star.bam_file
    File bam_index=star.bam_index
    File transcriptome_bam=star.transcriptome_bam
    File chimeric_junctions=star.chimeric_junctions
    File chimeric_bam_file=star.chimeric_bam_file
    File read_counts=star.read_counts
    File junctions=star.junctions
    File junctions_pass1=star.junctions_pass1
    Array[File] logs=star.logs

    #rnaseqc
    File gene_tpm=rnaseqc2.gene_tpm
    File gene_counts=rnaseqc2.gene_counts
    File exon_counts=rnaseqc2.exon_counts
    File metrics=rnaseqc2.metrics
    File insertsize_distr=rnaseqc2.insertsize_distr

    #rsem
    File genes=rsem.genes
    File isoforms=rsem.isoforms

  }
}

