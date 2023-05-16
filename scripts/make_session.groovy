import org.yaml.snakeyaml.*
import java.nio.file.*
import groovy.xml.*
import gngs.*

Utils.configureSimpleLogging()

log = java.util.logging.Logger.getLogger('make_session')

cli = new Cli(usage: 'make_session.groovy -vcf <vcf> -bam <bam> [-bam <bam>]... <session file name>')
cli.with {
    id 'id used to name the session', args:1, type: String, required: true
    vcf 'VCF file to load in session', args:1, type: File, required: true
    abs 'write absolute file paths'
    map 'Mapping of file prefix to other file prefix', args: '*', required: false
    bam 'RNA BAM file(s) to load in session', args:'*', type: File, required: true
}


opts = cli.parse(args)

def makeSession(String fam, File vcf, Map<String, File> bams) {
    def writer = new StringWriter()
    def xml = new MarkupBuilder(writer)
    def total_height = 0;
    def heights = [];
    xml.Session(genome: "hg38", hasGeneTrack: "true", hasSequenceTrack: "true", version: "8") {
        Resources {
            Resource(index:"https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz.tbi", 
                     name:"Gene",
                     path:"https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz",
                     type:"refgene")
            Resource(name: vcf, path: vcf)
            bams.each { k, v ->
                Resource(name: v, path: v)
            }
        }
        Panel(height: "120", name: "DataPanel") {
            Track(attributeKey: vcf, clazz: "org.broad.igv.variant.VariantTrack",
                  id: vcf, name: vcf, visible: "true")
            total_height += 120;
            heights << 120;
        }
        bams.each { nm, bam ->
            Panel(height: "300", name: nm) {
                Track(attributeKey: nm + " Coverage", autoScale: "true",
                      clazz: "org.broad.igv.sam.SpliceJunctionTrack",
                      id: bam + "_coverage", name: nm + " Coverage",
                      snpThreshold: "0.2", height: "30", visible: "true")
                Track(attributeKey: nm + " Junctions", autoScale: "true",
                      clazz: "org.broad.igv.sam.CoverageTrack",
                      id: bam + "_junctions", name: nm + " Junctions",
                      height: "60", visible: "true")
                Track(attributeKey: bam + "", autoScale: "true",
                      clazz: "org.broad.igv.sam.AlignmentTrack",
                      id: bam, name: bam, height: "220",
                      displayMode: "SQUISHED", visible: "true") {
                    RenderOptions()
                }
            }
            total_height += 300;
            heights << 300;
        }
        Panel(height: "60", name: "FeaturePanel") {
            Track(attributeKey: "Reference sequence",
                  clazz: "org.broad.igv.track.SequenceTrack",
                  id: "Reference sequence", name: "Reference sequence",
                  sequenceTranslationStrandValue: "POSITIVE",
                  shouldShowTranslation: "false", visible: "true")
            Track(attributeKey: "Gene",
                  clazz: "org.broad.igv.track.FeatureTrack",
                  id: "hg38_genes", name: "Gene",
                  height: "35", visible: "true")
        }
        total_height += 60;
        heights << 60;
        def running_total = 0
        def fracs = [];
        for (def i = 0; i < heights.size() - 1; i++) {
            running_total += heights[i];
            fracs << (running_total * 1.0 / total_height).toString();
        }
        PanelLayout(dividerFractions: fracs.join(','));
        HiddenAttributes {
            Attribute(name: "DATA FILE")
            Attribute(name: "DATA TYPE")
            Attribute(name: "NAME")
        }
    }
    return writer.toString();
}

if(!opts.arguments().size()) {
   cli.usage()
   throw new IllegalArgumentException("Please provide output file as argument")
}

def vcf = opts.vcf
def bams = opts.bams

if(opts.abs)  {
   vcf = vcf.absoluteFile
   bams = bams*.absoluteFile
}

bamMap = bams.collectEntries { bamFile ->
    [ new SAM(bamFile.path).samples[0], bamFile.path ]
}

log.info "BAM files / samples are: $bamMap"

if(opts.maps) {
   for(map in opts.maps) {
       def parts = map.tokenize('=')
       if(parts.size()!=2)
           throw new IllegalArgumentException("Please format map argument as <from>=<to>")
       vcf = new File(vcf.path.replaceAll(parts[0],parts[1]))
       bamMap = bamMap.collectEntries { sample, path ->
            [sample, path.replaceAll(parts[0],parts[1])]
       }
   }
}

result = makeSession(opts.id, vcf, bamMap)

outputFile = new File(opts.arguments()[0])
outputFile.text = result + '\n'

log.info "Wrote $outputFile.absolutePath"

