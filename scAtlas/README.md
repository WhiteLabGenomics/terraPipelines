The optimus pipeline supports the processing of any 3' single-cell and single-nucleus count data generated with the 10x Genomics v2 or v3 assay 
and can be applied to human (GRCh38) or mouse (GRCm38.p6) single-cell or single-nucleus datasets.

The smartseq2_single_sample pipeline processes a single sample (cell) and can be applied to human or mouse, stranded or unstranded, PE or SE, and 
plate- or fluidigm-based Smart-sew2 data.

The smartseq2 multi-sample pipeline is a wrapper around the standard single-sample pipeline that processes multiple cells by importing and running 
the Smart-seq2 Single Sample workflow for each cell (sample) and then merging the resulting Loom matrix output into a single Loom matrix 
containing raw counts and TPMs.

Source: https://broadinstitute.github.io/warp/

# test Optimus pipeline
cromwell run Optimus.wdl --inputs test_inputs/inputs_8k_pbmc.json

