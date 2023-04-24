version 1.0

import "fastqc.wdl" as fastqc_v1
import "multiqc.wdl" as multiqc_v1
import "star.wdl" as star_v1
import "rnaseqc2.wdl" as rnaseqc2_v1
import "rsem.wdl" as rsem_v1


workflow RNA_preprocessing_pipeline {

  input {

    File fastq1
    File fastq2
    String sample_id

    #rannotation GTF
    File genes_gtf="gs://ccle_default_params/references_gtex_gencode.v29.GRCh38.ERCC.genes.collapsed_only.gtf"

    #star_v1 index
    File star_index="gs://ccle_default_params/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v29_oh100.tar.gz"

    #rsem index
    File rsem_reference="gs://ccle_default_params/rsem_reference_GRCh38_gencode29_ercc.tar.gz"

  }

  call fastqc_v1.fastqc as fastqc1 {
    input:
      seqFile=fastq1,
      outdirPath=".",
  }

  call fastqc_v1.fastqc as fastqc2 {
    input:
      seqFile=fastq2,
      outdirPath=".",
  }

  call star_v1.star as star {
    input:
      prefix=sample_id,
      fastq1=fastq1,
      fastq2=fastq2,
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
      is_stranded="false",
      paired_end="true"
  }


  output {
    #fastqc
    File htmlReport1 = fastqc1.htmlReport
    File reportZip1 = fastqc1.reportZip
    File htmlReport2 = fastqc2.htmlReport
    File reportZip2 = fastqc2.reportZip
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

