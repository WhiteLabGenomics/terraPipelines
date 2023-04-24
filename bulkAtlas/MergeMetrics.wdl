task MergeMetrics {
  input {
    File alignment_summary_metrics
    File insert_size_metrics
    File picard_rna_metrics
    File duplicate_metrics
    File rnaseqc2_metrics
    File? fingerprint_summary_metrics
    String output_basename

    String docker =  "python:3.8-slim"
    Int cpu = 1
    Int memory_mb = 3000
    Int disk_size_gb = 10
  }

  String out_filename = output_basename + ".unified_metrics.txt"

  command <<<

    #
    # Script transpose a two line TSV
    #
    cat <<-'EOF' > transpose.py
    import csv, sys

    rows = list(csv.reader(sys.stdin, delimiter='\t'))

    for col in range(0, len(rows[0])):
      key = rows[0][col].lower()
      value = rows[1][col]
      if value == "?":
        value = "NaN"
      if key in ["median_insert_size", "median_absolute_deviation", "median_read_length", "hq_median_mismatches"]:
        value = str(int(float(value)))
      print(f"{key}\t{value}")
    EOF

    #
    # Script clean the keys, replacing space, dash and forward-slash with underscores,
    # and removing comma, single quote and periods
    #
    cat <<-'EOF' > clean.py
    import sys

    for line in sys.stdin:
      (k,v) = line.strip().lower().split("\t")
      transtable = k.maketrans({' ':'_', '-':'_', '/':'_', ',':None, '\'':None, '.' : None})
      print(f"{k.translate(transtable)}\t{v}")
    EOF

    # Process each metric file, transposing and cleaning if necessary, and pre-pending a source to the metric name

    echo "Processing Alignment Summary Metrics - Only PAIR line"
    cat ~{alignment_summary_metrics} | egrep "(CATEGORY|^PAIR)" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "picard_" $0}' >> ~{out_filename}

    echo "Processing Insert Size Metrics - removing various WIDTH metrics"
    cat ~{insert_size_metrics} | grep -A 1 "MEDIAN_INSERT_SIZE" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP|WIDTH)" | awk '{print "picard_" $0}' >> ~{out_filename}

    echo "Processing Picard RNA Metrics"
    cat ~{picard_rna_metrics} | grep -A 1 "RIBOSOMAL_BASES" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "picard_rna_metrics_" $0}' >> ~{out_filename}

    echo "Processing Duplicate Metrics"
    cat ~{duplicate_metrics} | grep -A 1 "READ_PAIR_DUPLICATES" | python transpose.py | awk '{print "picard_" $0}' >> ~{out_filename}

    echo "Processing RNASeQC2 Metrics"
    cat ~{rnaseqc2_metrics} | python clean.py | awk '{print "rnaseqc2_" $0}' >> ~{out_filename}

    if [[ -f "~{fingerprint_summary_metrics}" ]];
    then
      echo "Processing Fingerprint Summary Metrics - only extracting LOD_EXPECTED_SAMPLE"
      cat ~{fingerprint_summary_metrics} | grep -A 1 "LOD_EXPECTED_SAMPLE" | python transpose.py | grep -i "LOD_EXPECTED_SAMPLE" | awk '{print "fp_"$0}' >> ~{out_filename}
    else
      echo "No Fingerprint Summary Metrics found."
      echo "fp_lod_expected_sample	" >> ~{out_filename}
    fi    >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File unified_metrics = out_filename
  }
}

