//vim: autoindent ts=4 sw=4
import org.yaml.snakeyaml.*
import gngs.*
import groovy.json.*
import org.gitlab4j.api.*
import org.gitlab4j.api.models.*

options {
    enable_ase 'Enable ASE analysis (default: false)'
}

batch=file('.').absoluteFile.name

header = { msg ->
    int width = 80
    println ""
    println "*" * width
    println "* " + msg.center(width-4) + " *"
    println "*" * width
    println ""
}

header "RNADX Pipeline for Batch $batch"

title "RNADX Pipeline for Batch $batch"

load 'config_rnadx_home.groovy'

// Load file with paths to necessary software
load 'software_paths.groovy'

// Load site configuration last, so that it can override
// definitions, if necessary.
//
load 'config.groovy'

// Load bpipe stages required for this pipeline
//
load 'rnaseq_stages.groovy'
load 'common_stages.groovy'
load 'metadata.groovy'
/*
requires GITLAB_TOKEN : 'This pipeline requires a gitlab token to be configured to read sample metadata from gitlab'

gitlab = new GitLabApi("https://git.mcri.edu.au", f.Gitlab_Token);

cgst = gitlab.projectApi.getProjects('clingen_sample_tracking')[0]

println "Resolved gitlab project as $cgst.id ($cgst.title)"

*/

ESSENTIAL_R = "$RNADX_BASE/singularity/dockanomics_essential_r_4_0.sif"
PORTCULLIS_SIF = "$RNADX_BASE/singularity/maplesond_portcullis_latest.sif"
SALMONINDEX = "$RNADX_BASE/reference/hg38/salmon/salmon.genc.v38.index"

GIP_HTTP = "https://vcgsdev.mcri.edu.au:8443"

var FRASER_CONTROLS : false,
    DE_CONTROLS : false,
    RUN_FASTQC : true


header "Loading metadata ..."

if(file(args[0]).exists()) {
    meta = readMetaFromFile(args[0])
}
else {
    meta = readMetaFromAPI(args)
}


println("meta = " + JsonOutput.prettyPrint(JsonOutput.toJson(meta)))

fastqs = meta.samples.collectEntries { [it.key, it.value.fastqs ] }

samples = meta.samples*.key

println("Using hg38 genome for analysis....") 

CWD = System.getProperty("user.dir")

do_in_silico_depletion = {
    output.dir = "depleted"

    from("*_R1.fastq.gz","*_R2.fastq.gz") {
        def keep_names = inputs.collect { it.name.replaceAll('.fastq.gz','.keep.fastq.gz') }
        def toss_names = inputs.collect { it.name.replaceAll('.fastq.gz','.toss.fastq.gz') }
        produce(keep_names + toss_names) {
            println "Depleting $inputs.gz -> $keep_names"
            exec """
                $RNADX_BASE/tools/bbmap/38.91/bbmap/bbduk.sh
                    in=$input1
                    in2=$input2
                    out=$output1
                    out2=$output2
                    outm=$output3
                    outm2=$output4
                    ref=$RNADX_BASE/reference/depletion/rRNAs.fasta
            ""","rstuff"
        }
        forward keep_names
    }
}

deplete_and_trim = {
    output.dir = "depleted"

    from("*_R1.fastq.gz","*_R2.fastq.gz") {
        def keep_names = [
          inputs.find { it.name.contains('_R1.fastq.gz') },
          inputs.find { it.name.contains('_R2.fastq.gz') }
        ]*.name*.replaceAll('.fastq.gz','.keep.fastq.gz')

        def toss_names = keep_names*.replaceAll('.keep.fastq.gz','.toss.fastq.gz')

        produce(keep_names + toss_names + [sample + '.trim.stats.txt', sample + '.depletion.stats.txt']) {
            println "Depleting $inputs.gz -> $keep_names"
            exec """

                set -o pipefail

                java -cp $RNADX_BASE/rnadx/tools/bpipe/0.9.11/lib/bpipe.jar:$RNADX_BASE/rnadx/pipeline/libs/groovy-ngs-utils.jar  
                    gngs.tools.SplitFASTQ
                    -s 0,1
                    -r1 ${inputs.grep { file(it).name.contains("_R1.")}.join(",")}
                    -r2 ${inputs.grep { file(it).name.contains("_R2.")}.join(",")} |
                $RNADX_BASE/tools/bbmap/38.91/bbmap/bbduk.sh
                    -Xmx1g
                    threads=2
                    stats=$output.trim.stats.txt
                    in=stdin.fq
                    int=t
                    ktrim=r k=23 mink=11 hdist=1 tpe tbo
                    ref=$RNADX_BASE/reference/adapters/adapters.fa
                    out=stdout.fq |
                $RNADX_BASE/tools/bbmap/38.91/bbmap/bbduk.sh
                    -Xmx20g
                    threads=${threads-2}
                    in=stdin.fq
                    int=t
                    stats=$output.depletion.stats.txt
                    out=$output1
                    out2=$output2
                    outm=$output3
                    outm2=$output4
                    ref=$RNADX_BASE/reference/depletion/rRNAs.fasta
            ""","bbduk"
        }
        meta.samples[sample].kept = keep_names.collect { "depleted/" + it }
        meta.samples[sample].tossed = toss_names.collect { "depleted/" + it }
        forward keep_names.collect { "depleted/" + it }
    }
}

produce_batch_report = {
    output.dir = "reports"

    produce("sample-bams.tsv", "batch-report.html") {
        new File('reports/sample-bams.tsv').text = inputs.bam.collect {
            [ new SAM(it.toString()).samples[0], it.toString() ].join('\t') + '\n'
        }.sum()

        exec """

            echo "Creating batch report from $input.counts.txt"

            singularity exec -C
                -B $RNADX_BASE:/rnadx
                -H $CWD:/batch
                $ESSENTIAL_R
                R -e 'rmarkdown::render("/rnadx/rnadx/templates/batch-report.Rmd", 
                                        "html_document", 
                                         params=list(sample_bams="${output['sample-bams.tsv']}"),
                                         output_file="$output.html", 
                                         output_dir="/batch/reports", knit_root_dir=getwd())'
        ""","rstuff"
    }
}

run_portcullis = {
    output.dir = "portcullis"

    produce(branch.sample) {
        exec """
            singularity exec -C
                -B $RNADX_BASE:/rnadx
                -H $CWD:/batch
                $PORTCULLIS_SIF
                portcullis full -o $output /rnadx/reference/hg38/genome/GRCh38.primary_assembly.genome.fa $input.bam
        ""","rstuff"
    }
}

set_sample = {
    branch.sample = branch.name
    forward(meta.samples[sample].files)
}

set_kept = {
    branch.sample = branch.name
    branch.output_prefix = "";
    forward(meta.samples[sample].kept)
}

set_kept_prefix = {
    branch.output_prefix = "";
}

set_tossed = {
    branch.sample = branch.name
    branch.output_prefix = "removed_";
    forward(meta.samples[sample].tossed)
}

set_tossed_prefix = {
    branch.output_prefix = "removed_";
}

store_sample_bam = {
    meta.samples[branch.sample].bam = input.bam
    forward inputs
}

create_fraser_sample_sheets = {

    output.dir = 'fraser'

    produce("${sample_group}.test.tsv","${sample_group}.control.tsv") {
        println "Creating bam sample sheet $output.tsv from $inputs.bam ..."
        if(!file(output.test.tsv).exists()) {
            file(output.test.tsv).text = "sample\tbam\n" + sample_ids.collect { sample_id ->
               [sample_id, meta.samples[sample_id].bam].join('\t') + '\n'
            }.sum()

            file(output.control.tsv).text = "sample\tbam\n" + control_groups[sample_group].collect { sample_id, bam ->
               [sample_id, bam].join('\t') + '\n'
            }.sum()
        }
    }
}

create_sample_report = {

   branch.sample = branch.name

   requires REFSEQ : 'RefGene file downlaoded from UCSC',
            GENOME_FASTA : 'Human genome reference sequence'

   output.dir = "reports/html/$sample"

   def proband = meta.samples*.value.find { meta ->
       (sample in meta.parents) ||  // this sample is one of the parents
       (meta.parents && meta.identifier == sample) // this sample has parent
   }?.identifier

   if(!proband)
       fail "Could not find proband for family of $sample"

   def input_files = ["fraser.significant.results.tsv", "${sample}.results.tsv","${proband}.igv.session.xml"]

   def rare_vcf = meta.samples[proband].rare_vcf
   if(rare_vcf) {
       input_files << rare_vcf
   }

   from(input_files) produce("index.html") {
       def vcf_flag = rare_vcf ? "-vcf $input.vcf" : ""
       exec """
           $GROOVY -cp $GNGS_JAR $RNADX_BASE/rnadx/tools/rnadxview/scripts/CreateReportAssets.groovy
                -session $input.xml
                -fraser $input1
                -deseq $input2
                -refgene $REFSEQ
                -sample $sample $vcf_flag
                -ref $GENOME_FASTA
                -dir $output.dir
                -dist $RNADX_BASE/rnadx/tools/rnadxview/dist
       """
   }
}

make_session_file = {

       output.dir = "reports"

       var familyId : branch.name

       def members = meta.samples*.value.findAll { 
           it.familyId == familyId
       }

       def vcf = members.find { it.full_vcf }?.full_vcf
       if(!vcf)
           fail "No VCF could be found for family $familyId"

       def bams = members*.bam

       def proband = members.find { it.parents }?.identifier
       if(!proband)
           return "Not generating viewable report for ${members*.identifier} because no proband could be identified"

       from([vcf, *bams]) produce(proband+'.igv.session.xml') {
           exec """
               $GROOVY -cp $GNGS_JAR $RNADX_BASE/rnadx/scripts/make_session.groovy
                    -id $proband
                    -vcf $input1 ${inputs.bam.withFlag('-bam')}
                    -map /misc/vcgs=/Volumes/vcgs
                    -abs
                    $output.session.xml
           """
       }
}


control_groups = meta.control_groups.collectEntries { name, sample_sheet ->
    def control_base = file(sample_sheet).parentFile.absolutePath
    return [
        name,
        new TSV(sample_sheet).toListMap()*.bam
            .collect {
                "$control_base/$it"
            }
            .collectEntries { [ new SAM(it).samples[0], it ] }
    ]
}

all_fastq_files = meta.samples*.value*.files.flatten()

all_families = meta.samples*.value*.familyId.unique()

println "Families are: " + all_families

header "Scanning Available Control Groups ..."

println "Control groups are: \n" + JsonOutput.prettyPrint(JsonOutput.toJson(control_groups))

header "Assigning samples to control groups"

sample_groups = meta.samples*.value.groupBy  { sampleInfo ->
   def age = sampleInfo.parents ? "child" : "adult"
   return age + '_' + sampleInfo.sex
}

println "Control groups assignments are: \n" + JsonOutput.prettyPrint(JsonOutput.toJson(sample_groups))

probands = meta.samples.grep {
    it.value.parents && !it.value.parents.isEmpty()
}*.value

println "\nThe probands are: \n" + probands*.identifier.join('\n')

create_report_index = {

    send report('templates/report-index.html') to file: "reports/html/${batch}.html"
}


init_sample_group = {

    branch.sample_group = sample_group

    branch.sample_ids = sample_groups[sample_group]*.identifier

    println "Samples for group $sample_group are " +  sample_ids

    def group_bams = [
       *sample_groups[sample_group]*.bam,
       *control_groups[sample_group]*.value
    ]

    println "BAMs for sample group $sample_group are:\n\n" + 
            group_bams.join("\n") + "\n"

    forward(group_bams)
}

init_ase_proband = {

   def proband = probands.find { it.identifier == branch.name }

   branch.ASE_PROBAND = proband.identifier + '-' + proband.familyId

   if(!proband.full_vcf) 
       succeed "Skipping ASE analyiss for $ASE_PROBAND because this sample has no full trio VCF available"

   println "Resolved proband BAM file $proband.bam and VCF $proband.full_vcf for ASE on $ASE_PROBAND"

   forward(proband.full_vcf, proband.bam)
}

ase_disabled = {
   println "ASE was not enabled for this pipeline"
}

if(opts.enable_ase) {
    header "Precomputing regions for analysis ...."
    exome_regions = bed(REFFLAT_REGIONS).split(10)
}

run(all_fastq_files) {
    samples * [ set_sample + [ fastqc.when { RUN_FASTQC }, deplete_and_trim ] ] +
    samples * [ set_kept + star_map_1pass_sample_PE + star_map_2pass_sample_PE + store_sample_bam + index_bam ] +
    set_kept_prefix + 
    [
        sample_groups*.key.collect { sample_group ->

            init_sample_group.using(sample_group:sample_group) + 
          
            create_fraser_sample_sheets + [

                    // Run Fraser
                     fraser_metal,

                    // Count reads for DESeq
                    count_reads_RNA.using(sample_set: sample_group) +

                    sample_groups[sample_group]*.identifier * [ set_sample + deseq_simple ], // + produce_batch_report,

                    //# Steps to process the reads removed by depletion:
                    // #samples * [ set_tossed + star_map_1pass_sample_PE + star_map_2pass_sample_PE + sort_bam + index_bam ] +
                    //#    set_tossed_prefix + count_reads_RNA
            ] 
        } + merge_fraser_results,

        opts.enable_ase ?  probands*.identifier * [
            init_ase_proband + extract_ase_proband + exome_regions * [ count_ase_reads ] >>>  merge_ase_counts
        ] : ase_disabled
    ] +
    all_families * [ make_session_file ] +
    samples * [ create_sample_report ] + create_report_index
}
