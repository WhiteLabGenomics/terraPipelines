version 1.0

import "tasks/SplitNCigarReads.wdl" as SplitNCigarReads_task
import "tasks/gtfToCallingIntervals.wdl" as gtfToCallingIntervals_task
import "tasks/MarkDuplicates.wdl" as MarkDuplicates_task
import "tasks/BaseRecalibrator.wdl" as BaseRecalibrator_task
import "tasks/ApplyBQSR.wdl" as ApplyBQSR_task
import "tasks/ScatterIntervalList.wdl" as ScatterIntervalList_task
import "tasks/HaplotypeCaller.wdl" as HaplotypeCaller_task
import "tasks/MergeVCFs.wdl" as MergeVCFs_task
import "tasks/VariantFiltration.wdl" as VariantFiltration_task
import "tasks/Funcotate.wdl" as Funcotate_task

## Copyright Broad Institute, 2019
##
## Workflows for processing RNA data for germline short variant discovery with GATK (v4) and related tools
##
## Requirements/expectations :
## - BAM 
##
## Output :
## - A BAM file and its index.
## - A VCF file and its index. 
## - A Filtered VCF file and its index. 
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs. 

workflow call_variants {

    input {
        File inputBam
        String sampleName

        File refFasta
        File refFastaIndex
        File refDict

        String gatk4_docker = "broadinstitute/gatk:4.4.0.0"
        String gatk_path = "/gatk/gatk"
        String star_docker = "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-40ead6e"

        Array[File] knownVcfs
        Array[File] knownVcfsIndices

        File? dbSnpVcf
        File? dbSnpVcfIndex

        Int? minConfidenceForVariantCalling

        ## Inputs for STAR
        Int? readLength
        File? zippedStarReferences
        File genes_gtf

        ## Optional user optimizations
        Int haplotypeScatterCount=6

        Int preemptible_tries=3
        # funcotator
        String? sequencing_center
        String? sequence_source
        String? funco_reference_version
        String funco_output_format="VCF"
        Boolean funco_compress=true
        Boolean funco_use_gnomad_AF=false
        File? funco_data_sources_tar_gz
        String? funco_transcript_selection_mode
        File? funco_transcript_selection_list
        Array[String]? funco_annotation_defaults
        Array[String]? funco_annotation_overrides
        Array[String]? funcotator_excluded_fields
        Boolean funco_filter_funcotations=true
        String? funcotator_extra_args
    }

    Int funco_tar_size = if defined(funco_data_sources_tar_gz) then ceil(size(funco_data_sources_tar_gz, "GB") * 3) else 100
    
    call gtfToCallingIntervals_task.gtfToCallingIntervals as gtfToCallingIntervals {
        input:
            gtf = genes_gtf,
            ref_dict = refDict,
            preemptible_count = preemptible_tries,
            gatk_path = gatk_path,
            docker = gatk4_docker
    }

    call MarkDuplicates_task.MarkDuplicates as MarkDuplicates {
        input:
            input_bam = inputBam,
            base_name = sampleName + ".dedupped",
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


    call SplitNCigarReads_task.SplitNCigarReads as SplitNCigarReads {
        input:
            input_bam = MarkDuplicates.output_bam,
            input_bam_index = MarkDuplicates.output_bam_index,
            base_name = sampleName + ".split",
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            ref_dict = refDict,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }

    # if dbSNP exists:
    Boolean canDoBQSR = defined(dbSnpVcf) && defined(dbSnpVcfIndex)
    if(canDoBQSR){
        call BaseRecalibrator_task.BaseRecalibrator as BaseRecalibrator {
            input:
                input_bam = SplitNCigarReads.output_bam,
                input_bam_index = SplitNCigarReads.output_bam_index,
                recal_output_file = sampleName + ".recal_data.csv",
                dbSNP_vcf = select_first([dbSnpVcf]),
                dbSNP_vcf_index = select_first([dbSnpVcfIndex]),
                known_indels_sites_VCFs = knownVcfs,
                known_indels_sites_indices = knownVcfsIndices,
                ref_dict = refDict,
                ref_fasta = refFasta,
                ref_fasta_index = refFastaIndex,
                preemptible_count = preemptible_tries,
                docker = gatk4_docker,
                gatk_path = gatk_path
        }

        call ApplyBQSR_task.ApplyBQSR as ApplyBQSR {
            input:
                input_bam =  SplitNCigarReads.output_bam,
                input_bam_index = SplitNCigarReads.output_bam_index,
                base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
                ref_fasta = refFasta,
                ref_fasta_index = refFastaIndex,
                ref_dict = refDict,
                recalibration_report = BaseRecalibrator.recalibration_report,
                preemptible_count = preemptible_tries,
                docker = gatk4_docker,
                gatk_path = gatk_path
        }
    }


    call ScatterIntervalList_task.ScatterIntervalList as ScatterIntervalList {
        input:
            interval_list = gtfToCallingIntervals.interval_list,
            scatter_count = haplotypeScatterCount,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


    scatter (interval in ScatterIntervalList.out) {
        call HaplotypeCaller_task.HaplotypeCaller as HaplotypeCaller {
            input:
                input_bam = select_first([ApplyBQSR.output_bam, SplitNCigarReads.output_bam]),
                input_bam_index = select_first([ApplyBQSR.output_bam_index, SplitNCigarReads.output_bam_index]),
                base_name = sampleName + ".hc",
                interval_list = interval,
                ref_fasta = refFasta,
                ref_fasta_index = refFastaIndex,
                ref_dict = refDict,
                dbSNP_vcf = dbSnpVcf,
                dbSNP_vcf_index = dbSnpVcfIndex,
                stand_call_conf = minConfidenceForVariantCalling,
                preemptible_count = preemptible_tries,
                docker = gatk4_docker,
                gatk_path = gatk_path
        }

    }

    call MergeVCFs_task.MergeVCFs as MergeVCFs {
        input:
            input_vcfs = HaplotypeCaller.output_vcf,
            input_vcfs_indexes =  HaplotypeCaller.output_vcf_index,
            output_vcf_name = sampleName + ".g.vcf.gz",
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }
    
    call VariantFiltration_task.VariantFiltration as VariantFiltration {
        input:
            input_vcf = MergeVCFs.output_vcf,
            input_vcf_index = MergeVCFs.output_vcf_index,
            base_name = sampleName + ".variant_filtered.vcf.gz",
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            ref_dict = refDict,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }

    Int disk_pad = 15
    Float file_size_multiplier=2.25

    call Funcotate_task.Funcotate as Funcotate {
        input:
            ref_fasta = refFasta,
            ref_fai = refFastaIndex,
            ref_dict = refDict,
            input_vcf = VariantFiltration.output_vcf,
            input_vcf_idx = VariantFiltration.output_vcf_index,
            # could be none
            reference_version = select_first([funco_reference_version, "hg38"]),
            output_file_base_name = basename(VariantFiltration.output_vcf, ".vcf") + ".annotated",
            output_format = funco_output_format,
            compress = funco_compress,
            use_gnomad = funco_use_gnomad_AF,
            data_sources_tar_gz = funco_data_sources_tar_gz,

            sequencing_center = sequencing_center,
            sequence_source = sequence_source,
            transcript_selection_mode = funco_transcript_selection_mode,
            transcript_selection_list = funco_transcript_selection_list,
            annotation_defaults = funco_annotation_defaults,
            annotation_overrides = funco_annotation_overrides,
            funcotator_excluded_fields = funcotator_excluded_fields,
            interval_list=gtfToCallingIntervals.interval_list, # maybe not the right one
            filter_funcotations = funco_filter_funcotations,
            extra_args = funcotator_extra_args,
            disk_space = ceil(size(VariantFiltration.output_vcf, "GB") * file_size_multiplier)  + funco_tar_size + disk_pad,
            gatk_docker=gatk4_docker,
            machine_mem=4000,
            preemptible=2,
            max_retries=2,
            cpu=2
    }

    output {
        File merged_vcf = MergeVCFs.output_vcf
        File merged_vcf_index = MergeVCFs.output_vcf_index
        File variant_filtered_vcf = VariantFiltration.output_vcf
        File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
        File funcotated_output_file = Funcotate.funcotated_output_file
        File funcotated_output_file_index = Funcotate.funcotated_output_file_index
    }
}