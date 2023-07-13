# Optimus
The optimus pipeline supports the processing of any 3' single-cell and single-nucleus count data generated with the 10x Genomics v2 or v3 assay and can be applied to human (GRCh38) or mouse (GRCm38.p6) single-cell or single-nucleus datasets.

# Smart-seq2
The smartseq2_single_sample pipeline processes a single sample (cell) and can be applied to human or mouse, stranded or unstranded, PE or SE, and plate- or fluidigm-based Smart-seq2 data.

The smartseq2 multi-sample pipeline is a wrapper around the standard single-sample pipeline that processes multiple cells by importing and running the Smart-seq2 Single Sample workflow for each cell (sample) and then merging the resulting Loom matrix output into a single Loom matrix containing raw counts and TPMs.

source: https://broadinstitute.github.io/warp/

## Test optimus pipeline:
# validate wdl file
womtool validate Optimus.wdl

# output inputs file
womtool inputs Optimus.wdl > test_data/optimus_inputs.json

# download 10X test data
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3

# copy test data to WLG GCP
# https://cg.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_fastqs.tar

# copy reference and resource files to WLG GCP
cat ref-files.txt | gsutil -m cp -I gs://whitelabgx-references/hg38/
cat resource-files.txt | gsutil -m cp -I gs://whitelabgx-references/resources/

# https://www.gencodegenes.org/human/release_27.html
gsutil cp gencode.v27.annotation.gtf gs://whitelabgx-references/hg38/

# test locally using cromwell
cromwell run Optimus.wdl --inputs test_data/8k_pbmc_v2_inputs.json 

## Test smartseq2 pipelines:
# validate wdl files and generate input files
womtool validate SmartSeq2SingleSample.wdl
womtool inputs SmartSeq2SingleSample.wdl > test_data/SmartSeq2SingleSample_inputs.json
womtool validate MultiSampleSmartSeq2.wdl
womtool inputs MultiSampleSmartSeq2.wdl > test_data/MultiSampleSmartSeq2_inputs.json