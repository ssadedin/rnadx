sra_to_fastq_PE = {
  // Convert SRA files to fastq format, paired-end
  output.dir = "fastq"
  transform(".sra") to ("_1.fastq.gz","_2.fastq.gz"){
    exec """
      $FASTQ_DUMP --split-3 --gzip -O $output.dir $input.sra;
    ""","small"
  }
}

sra_to_fastq_SE = {
  // Convert SRA files to fastq format, single-end
  output.dir = "fastq"
  transform(".sra") to (".fastq.gz"){
    exec """
      $FASTQ_DUMP --split-3 --gzip -O $output.dir $input.sra;
    ""","small"
  }
}

concat_fastq_SE = {
    // Concatenate fastq files from the same sample
    doc "Concatenate fastq files from the same run"
    output.dir = "catfastq"

    produce(output.gz.prefix.prefix.prefix + ".cat.fastq.gz"){
        exec """
                cat $inputs.gz > $output.gz
        ""","tiny"
    }
}

concat_fastq = {
    // Concatenate fastq files from the same sample
    doc "Concatenate fastq files from the same run"
    output.dir = "catfastq"

    produce(output.gz.prefix.prefix.prefix + ".fastq.gz"){
        exec """
                cat $inputs.gz > $output.gz
        ""","tiny"
    }
}

fastqc = {
    // fastqc quality control
    doc "Quality control using FASTQC"
    output.dir = "fastqc"

   from('*.fastq.gz') {

       def output_names = inputs.collect { it.name.replaceAll('.fastq.gz','_fastqc.zip') }

       produce(output_names) {
            exec """
            $FASTQC -t $threads -o $output.dir $inputs.gz
        ""","medium"
       }
    }
    //println inputs;
    forward inputs
}

multiqc = {
   // summarise statistics from all tools using multiqc
   doc "Pipeline summary report by MultiQC"
   output.dir = "mulitQC"
   exec """
	multiqc -o $output.dir . --ignore .bpipe
   ""","multiqc"
}

trim_SE = {
    // trim single-end reads using trimmomatic
    doc "Trim poor quility bases and/or adapter sequences from reads using Trimmomatic"
    output.dir = "trimmed"

    filter("trim"){
        exec """
            $TRIMOMMATIC SE -threads $threads $input.gz $output.gz $CLIPSTRING
            TRAILING:$TRAILING LEADING:$LEADING MINLEN:$MINLEN
        ""","trimmomatic"
    }
}

trim_PE = {
   // trim paired-end reads using trimmomatic
   doc "Trim poor quility bases and/or adapter sequences from reads using Trimmomatic"
   output.dir="trimmed"

   from('*.fastq.gz') {

       def output_names = inputs.collect { it.name.replaceAll('.fastq.gz','.trim.fastq.gz') }

       produce(output_names) {
         exec """
                $TRIMOMMATIC PE -threads $threads
                $input1.gz $input2.gz
                $output1.gz ${output1.prefix}.unpaired.gz
                $output2.gz ${output2.prefix}.unpaired.gz ILLUMINACLIP:$adapters:$ILLUMINACLIP
                LEADING:$LEADING TRAILING:$TRAILING MINLEN:$MINLEN
         ""","trimmomatic"

         println "The trimmed output for $sample is: $output1.gz, $output2.gz"

        // Is this just voodoo?
        def sample = "$sample"
        // def sample = file(input1.gz).name.toString().tokenize('.')[0]

         meta[sample].files.addAll(0,[output1.gz, output2.gz])
     }
   }
}

sort_bam = {
    // sort bam files
    doc "Sort bam files"
    output.dir = "${branch.output_prefix}sorted"

    filter("srt") {
        exec """
		$SAMTOOLS sort -m 3G -@ $threads -o $output $input
	""","srtindex"
    }
}

index_bam = {
    doc "Index a bam file"
    output.dir = "${branch.output_prefix}sorted"

    transform("bam") to ("bam.bai") {
        exec """
          $SAMTOOLS index $input.bam
        """
    }
    forward input
}

count_mapped_SE = {
    // Count total number of single-end mapped reads in bam file (no multi-map)
    doc "Count total number of single-end mapped reads in bam file (no multi-map)"
    output.dir = "${branch.output_prefix}sorted"

    transform("txt"){
        exec """
                $SAMTOOLS view -F 0x904 -c $input > $output
        ""","tiny"
    }
}

count_mapped_PE = {
    // Count total number of mapped paired-end reads in bam file (no multi-map)
    doc "Count total number of mapped paired-end reads in bam file (no multi-map)"
    output.dir = "${branch.output_prefix}stats"

    transform("txt"){
        exec """
		$SAMTOOLS view -F 0x4 $input.bam | cut -f 1 | sort | uniq | wc -l > $output
        ""","tiny"
    }
}

coverage_by_region = {
    // Get read coverage for each base across a particular region(s)
    doc "Get read coverage for each base across a particular region(s)"
    output.dir = "${branch.output_prefix}stats"

    transform("bed"){
        exec """
                $BEDTOOLS coverage -d -split -abam $input.bam -b $REGIONS_BED > $output.bed
        ""","regioncov"
    }
}

hard_link_vcf = {
    output.dir = "vcfs"
    produce(file(input).name){
        exec """
            ln $input $output
        ""","trivial"
    }
}
