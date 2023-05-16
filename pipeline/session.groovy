import org.yaml.snakeyaml.*
import java.nio.file.*
import groovy.xml.*

// Load file with paths to necessary software
load 'software_paths.groovy'

// Load bpipe stages required for this pipeline
//
load 'rnaseq_stages.groovy'
load 'common_stages.groovy'

// Load site configuration last, so that it can override
// definitions, if necessary.
//
load 'config.groovy'

templateFilename = "$RNADX_BASE/rnadx/scripts/igv-session-template.xml".toString()

familyMeta = file(args[0]).withReader { r -> Collections.synchronizedMap(new Yaml().load(r)) }

families = familyMeta*.key

set_rare_vcf = {
    def fam = branch.name
    def vcfPath = familyMeta[fam]['vcf-rare']
    forward(vcfPath)
}

def makeSession(fam, vcf, bams) {
    def writer = new StringWriter()
    def xml = new MarkupBuilder(writer)
    def total_height = 0;
    def heights = [];
    xml.Session(genome: "hg38", hasGeneTrack: "true", hasSequenceTrack: "true", version: "8") {
        Resources {
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
                Track(attributeKey: bam + " Coverage", autoScale: "true",
                      clazz: "org.broad.igv.sam.SpliceJunctionTrack",
                      id: bam + "_coverage", name: bam + " Coverage",
                      snpThreshold: "0.2", height: "30", visible: "true")
                Track(attributeKey: bam + " Junctions", autoScale: "true",
                      clazz: "org.broad.igv.sam.CoverageTrack",
                      id: bam + "_junctions", name: bam + " Junctions",
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
        Panel(height: "40", name: "FeaturePanel") {
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
        total_height += 40;
        heights << 40;
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

make_session_file = {
    def fam = branch.name

    println "making session file for ${fam}"
    produce("${fam}_igv_session.xml") {
        def vcf = "${input}".toString()
        def members = familyMeta[fam]['members']
        def bams = familyMeta[fam]['rna-bams']
        def dict = GroovyCollections.transpose([members, bams]).collect { pair -> [(pair[0]): pair[1]]}.collectEntries { it }
        def f = file("${fam}_igv_session.xml");
        f.text = makeSession(fam, vcf, dict)
    }
}

run() {
    families * [ set_rare_vcf + hard_link_vcf + make_session_file ]
}

