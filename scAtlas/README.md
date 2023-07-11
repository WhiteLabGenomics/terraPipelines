# Optimus
The optimus pipeline supports the processing of any 3' single-cell and single-nucleus count data generated with the 10x Genomics v2 or v3 assay 
and can be applied to human (GRCh38) or mouse (GRCm38.p6) single-cell or single-nucleus datasets.

# Smart-seq2
The smartseq2_single_sample pipeline processes a single sample (cell) and can be applied to human or mouse, stranded or unstranded, PE or SE, and 
plate- or fluidigm-based Smart-seq2 data.

The smartseq2 multi-sample pipeline is a wrapper around the standard single-sample pipeline that processes multiple cells by importing and running the Smart-seq2 Single Sample workflow for each cell (sample) and then merging the resulting Loom matrix output into a single Loom matrix 
containing raw counts and TPMs.

Source: https://broadinstitute.github.io/warp/

# Test optimus pipeline:
# validate wdl file
womtool validate Optimus.wdl

# output inputs file
womtool inputs Optimus.wdl > test_inputs/optimus_inputs.json

# download 10X test data
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3

# test locally using cromwell
cromwell run Optimus.wdl --inputs test_inputs/10k_pbmc_v3.json
