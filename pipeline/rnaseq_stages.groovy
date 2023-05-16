star_genome_gen = {
    //Generate STAR genome index
    doc "Generate STAR genome index"

    produce(STAR_GENOME_DIR + "Genome") {
        exec """

            mkdir -p $STAR_GENOME_DIR

            $STAR --runMode genomeGenerate 
                --runThreadN $threads
                --genomeDir $STAR_GENOME_DIR
                --genomeFastaFiles $GENOME_FASTA
                --sjdbGTFfile $GENOME_ANNOTATION
                --sjdbOverhang $SJBOHANG
        ""","stargen"
    }
}

star_map_1pass_global_PE = {
    doc "Map paired-end reads using the STAR aligner: 1st pass"
    output.dir = "${branch.output_prefix}mapped"

    produce ("SJ.out.tab") {
        exec """
            $STAR --genomeDir $STAR_GENOME_DIR 
                --readFilesIn ${inputs.gz.join(',')}
                --readFilesCommand zcat 
                --outSAMtype None 
                --limitOutSJcollapsed 5000000
                --runThreadN $threads 
                --sjdbGTFfile $GENOME_ANNOTATION 
                --outFileNamePrefix ${output.dir}/
        ""","star1pass"
    }
}

star_map_2pass_global_PE = {
    //Map paired-end reads using the STAR aligner: 2nd pass
    doc "Map paired-end reads using the STAR aligner: 2nd pass"
    output.dir = "${branch.output_prefix}mapped"

    from("${sample}_SJ.out.tab") produce("${sample}.Aligned.sortedByCoord.out.bam") {
        exec """
            $STAR --genomeDir $STAR_GENOME_DIR 
                --readFilesIn $input1.gz $input2.gz
                --sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab 
                --sjdbOverhang $SJBOHANG
                --outSAMunmapped Within
                --outSAMattrRGline ID:$sample PL:Illumina LB:$sample SM:$sample
                --outFileNamePrefix ${output.prefix.prefix.prefix}. 
                --readFilesCommand zcat
                --limitSjdbInsertNsj 2000000
                --outSAMtype BAM SortedByCoordinate
                --runThreadN $threads
        ""","star2pass"

        println("The sample is $sample")
    }
}

star_map_1pass_sample_PE = {
    doc "Map paired-end reads using the STAR aligner: 1st pass"
    output.dir = "${branch.output_prefix}mapped/1pass"

    produce("${sample}_SJ.out.tab") {
        exec """
            $STAR --genomeDir $STAR_GENOME_DIR 
                --readFilesIn ${inputs.gz.join(',')}
                --readFilesCommand zcat 
                --outSAMtype None 
                --limitOutSJcollapsed 5000000
                --runThreadN $threads 
                --sjdbGTFfile $GENOME_ANNOTATION 
                --outFileNamePrefix ${output.dir}/${sample}_
        ""","star1pass"
    }
    forward inputs
}

star_map_2pass_sample_PE = {
    //Map paired-end reads using the STAR aligner: 2nd pass
    doc "Map paired-end reads using the STAR aligner: 2nd pass"
    output.dir = "${branch.output_prefix}sorted"

    from("${sample}_SJ.out.tab") produce("${sample}.Aligned.sortedByCoord.out.bam") {
        exec """
            $STAR --genomeDir $STAR_GENOME_DIR 
                --readFilesIn $input1.gz $input2.gz
                --sjdbFileChrStartEnd $input.out.tab 
                --sjdbOverhang $SJBOHANG
                --outSAMunmapped Within
                --outSAMattrRGline ID:$sample PL:Illumina LB:$sample SM:$sample
                --outFileNamePrefix ${output.prefix.prefix.prefix.prefix}. 
                --readFilesCommand zcat
                --limitSjdbInsertNsj 2000000
                --outSAMtype BAM SortedByCoordinate
                --outSAMmultNmax 4
                --runThreadN $threads
        ""","star2pass"

        println("The sample is $sample")
    }
}

count_reads_RNA = {
    //Count reads across features using RNA data
    doc "Count reads across features from RNA data"

    var sample_set : 'test'

    output.dir = "${branch.output_prefix}counts"

    produce("${sample_set}.feature.counts.txt") {
        exec """
            $FEATURECOUNTS --primary -p -t exon -g gene_name -T $threads  
                -a $GENOME_ANNOTATION
                -s 0
                -o $output $inputs.bam
        ""","count"
    }
}

count_multi_reads_RNA = {
    //Count reads across features using RNA data
    doc "Count reads across features from RNA data including multimappers"
    output.dir = "${branch.output_prefix}counts"

    produce("feature.counts.multi.txt") {
    exec """
        $FEATURECOUNTS --primary -M -p -t exon -g gene_name -T $threads  
            -a $GENOME_ANNOTATION
            -s 0
            -o $output $inputs.bam
    ""","count"
    }
}

salmon_quant_PE = {
    doc "Quasi-map and quantify paired-end RNAseq reads using Salmon"
    output.dir = "${branch.output_prefix}salmon"

    exec """
        module load salmon &&
        $SALMON quant -i $SALMONINDEX -l A \
            -1 $input1.gz \
            -2 $input2.gz \
            -p $threads \
            --seqBias \
            --gcBias \
            --validateMappings \
            -o ${output.dir}/${sample}_quant/
    ""","salmon"
}

fraser_metal = {

    output.dir = "fraser"

    produce("fraser.${sample_group}.significant.results.tsv") {
        exec """

            export R_LIBS_USER=$RNADX_BASE/tools/r_libs

            export GENCODE_ANNOTATION_GTF=$RNADX_BASE/reference/hg38/annotation/gencode.v37.annotation.gtf

            export RNADX_BATCH=$batch

            export FRASER_CONTROLS="$input.control.tsv"

            export FRASER_TEST="$input.test.tsv"

            export FRASER_OUTPUT_FILE="$output.tsv"

            export FRASER_WORKING_DIR="$sample_group"

            export BATCH_BASE="$RNADX_BASE/batches/$batch"

            Rscript --vanilla $RNADX_BASE/rnadx/scripts/run-fraser.R

        ""","fraser_metal"
    }
}

fraser = {
    output.dir = "reports"

    produce("fraser.results.tsv") {
        exec """
            RNADX_BATCH=\$(basename $CWD) &&
            singularity exec -C
                --env GENCODE_ANNOTATION_GTF=/rnadx/reference/hg38/annotation/gencode.v37.annotation.gtf
                --env RNADX_BATCH=$batch
                --env BATCH_BASE="/batches/$batch"
                -B $RNADX_BASE:/rnadx
                -H $RNADX_BASE/batches:/batches
                $ESSENTIAL_R
                Rscript /rnadx/rnadx/scripts/run-fraser.R
        ""","fraser"
    }
}

merge_fraser_results = {

    output.dir = "fraser"

    produce("fraser.significant.results.tsv") {
        exec """
            java -cp $GROOVY_JAR:$GNGS_JAR gngs.tools.HeaderCat
            -o $output.tsv
            $inputs.significant.results.tsv
        ""","tiny"
    }
}

deseq_simple = {

    requires DESIMPLE_HOME : 'Location of the DESIMPLE repository'

    output.dir = "results"
    if (branch.output_prefix) {
        output.dir = output_prefix + output.dir
    }

    def f = file("${sample}.de-samplesheet.tsv")

    if(!f.exists()) {
        def ss = 'samplename\tfilename\n' + 
            inputs.bam.collect { bam ->
                def sample = new gngs.SAM(file(bam).path).samples[0]
                [sample, bam].join('\t')
            }.join('\n')   + '\n'

        f.text = ss + '\n'
    }

    produce("${sample}.results.tsv", "${sample}.results.pdf") {
        exec """
            $DESIMPLE_HOME/deseq2-simple 
               -S ${sample}.de-samplesheet.tsv
               -n  -s $sample
               -P $output.pdf
               -o $output.tsv $inputs.counts.txt
        ""","deseq"
    }
}

/**
 * Allele Specific Expression
 */

extract_ase_proband = { 

   requires REFFLAT_REGIONS : 'BED file containing transcripts for evaluation of ASE',
            ASE_PROBAND : 'The name of the target sample to analyse within the VCF provided'

   output.dir = "ase/$ASE_PROBAND"

   transform('vcf.gz') to("${ASE_PROBAND}.singleton.vcf.gz") {
       exec """
            set -o pipefail

            unset GROOVY_HOME ; export GROOVY="$GROOVY"

            $BCFTOOLS view -T $REFFLAT_REGIONS -s $ASE_PROBAND $input.vcf.gz | 
            $GNGS 'prev = null; VCF.filter { v -> 
               if(Region.isMinorContig(v.chr) || v.alts.contains("*") || (v.type != "SNP") || !v.het || prev?.pos==v.pos) return false; prev = v; }
            ' |
            $BGZIP -c > $output.vcf.gz

            $TABIX -p vcf $output.vcf.gz

        ""","extractvariants"
    }
}

count_ase_reads = {

    output.dir = "ase/$ASE_PROBAND"

    exec """

        java -Xmx${memory}g -jar $GATK SplitNCigarReads 
                         -R $GENOME_FASTA 
                         -I $input.bam 
                         -O $output.bam
                         -L $region.bed

        $SAMTOOLS index $output.bam

        java -Xmx${memory}g -jar $GATK ASEReadCounter 
                         -L $region.bed
                         -R ${GENOME_FASTA} 
                         -I $output.bam 
                         -V $input.vcf.gz
                         -O $output.tsv 
    ""","gatk_ase"
}

merge_ase_counts = {

    output.dir = "ase/$ASE_PROBAND"

    exec """
        unset GROOVY_HOME ; export GROOVY="$GROOVY"

        head -n 1 $input.tsv > $output.tsv

        ${GNGS}tool HeaderCat -c NONE -o - $inputs.tsv | tail -n +2 | sort -k1,1V -k2,2n -k3,3n >> $output.tsv
    """
}



