## RNAseq quantification using kallisto
## author: Dina Zielinski

## TASKS

task align_pe {
	File fastq_1
	File fastq_2
	File kallisto
	File index
	String outdir
	Int? frag_len
	Float? sd_frag_len
		
	command {
		${kallisto} quant \
		-i ${index} \
		-o ${outdir} \
		--plaintext \
		${fastq_1} ${fastq_2}
	}

	output {
		File quants = "${outdir}/abundance.tsv"
  		#Array[File] quants = glob("*_abundance.tsv")
	}

	#runtime { docker: 'ubuntu:latest' }
}

task align_se {
	File fastq_1
	File kallisto
	File index
	String outdir
	Int frag_len
	Float sd_frag_len
		
	## to do: write sep alignment script & call in command
	command {
		${kallisto} quant \
		-i ${index} \
		-o ${outdir} \
		--plaintext \
		--single \
		-l ${frag_len}
		-s ${sd_frag_len} \
		$fastq_1
	}

	output {
		File quants = "${outdir}/abundance.tsv"
	  	#Array[File] quants = glob("*_abundance.tsv")
	}
	
	#runtime { docker: 'ubuntu:latest' }

}

workflow rna {
	
	File inputSamplesFile
	Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)
	String outDir
	
	File kallisto
	File index
	Boolean PE
	Int? frag_len
	Float? sd_frag_len

	if (PE) { 	
		scatter (sample in inputSamples) {
			call align_pe { input:
				fastq_1 = sample[0],
				fastq_2 = sample[1],
				kallisto = kallisto,
				index = index,
				outdir = outDir,
				frag_len = frag_len,
				sd_frag_len = sd_frag_len
			}
		}
	}

	if (!PE) { 
		scatter (sample in inputSamples) {
			call align_se { input:
				fastq_1 = sample[0],
				kallisto = kallisto,
				index = index,
				outdir = outDir,
				frag_len = frag_len,
				sd_frag_len = sd_frag_len
			}
		}
	}
}
