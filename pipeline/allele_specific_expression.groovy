import org.yaml.snakeyaml.*
import java.nio.file.*
import groovy.xml.*
import gngs.*

load "config_rnadx_home.groovy"
load "metadata.groovy"

meta = readMetaFromFile(args[0])

// Load file with paths to necessary software
load 'software_paths.groovy'

// Load bpipe stages required for this pipeline
//
load 'common_stages.groovy'
load 'rnaseq_stages.groovy'

// Load site configuration last, so that it can override
// definitions, if necessary.
load 'config.groovy'

requires REFFLAT_REGIONS : 'File defining regions of Refseq coding regions'

sample_ids = channel(meta.samples*.value.collect { it.identifier + '-' + it.familyId}).named('ASE_PROBAND')

println "Processing samples: \n" + sample_ids

exome_regions = bed(REFFLAT_REGIONS).split(8)

set_ase_files = {
    def plainId = ASE_PROBAND.tokenize('-')[0]
    def familyId = meta.samples[plainId].familyId
    def vcf = meta.samples*.value
                  .grep { it.familyId == familyId }
                  .find { it.full_vcf }
                  .full_vcf

    println "Located vcf $vcf for sample $ASE_PROBAND in family $familyId"

    // For now we just try to find the BAM file based on the sample id inside it
    def bam = inputs.bam.find { 
      def check_bam = new SAM(it.toString())
      println "Checking bam $it : sample = " + check_bam.samples[0]
      return check_bam.samples[0] == plainId 
    }
    assert bam != null : "Unable to find BAM file for sample $ASE_PROBAND : was this provided in the inputs?"

    forward(vcf, bam)
}

run(meta.samples*.value*.full_vcf.grep { it } + args.grep { file(it).name.endsWith('.bam') }) {
    sample_ids * [
        set_ase_files + extract_ase_proband +
        exome_regions * [ count_ase_reads ] >>>  merge_ase_counts
    ] 
}
