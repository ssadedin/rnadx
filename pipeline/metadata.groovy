import org.yaml.snakeyaml.*
import gngs.*
import org.gitlab4j.api.*
import org.gitlab4j.api.models.*


var CONTROL_SET : false

Map<String, Map> readMetaFromFile(metaFile) {
    // treat the argument as a literal samples.yaml file to use
    def allMeta = file(metaFile)
            .withReader { r -> Collections.synchronizedMap(new Yaml().load(r)) };
    allMeta.families = [:]
    allMeta.rareVcfs = [:]
    allMeta.fullVcfs = [:]
    allMeta.samples.each {
        if (!allMeta.families.containsKey(it.familyId)) {
            allMeta.families[it.familyId] = []
        }
        allMeta.families[it.familyId] << it.identifier;
        if (it.rare_vcf) {
            allMeta.rareVcfs[it.familyId] = it.rare_vcf;
        }
        if (it.full_vcf) {
            allMeta.fullVcfs[it.familyId] = it.full_vcf;
        }
    }

    // Convert the samples property to be a Map keyed on the
    // the sample id for easy lookup
    allMeta.samples = allMeta.samples
                             .collectEntries { meta -> [
                                    meta.identifier, [
                                        *:meta,
                                        files : [*meta.fastqs],
                                    ]
                             ]};
    return allMeta;
}

def readMetaFromAPI(List samples) {

    if(file('samples.yaml').exists()) {
        println "\nNOTE: existing samples.yaml file found. Using these results for sample meta data.\n"
        return readMetaFromFile('samples.yaml')
    }

    // treat the argument as a sample id to analyse
    def ws = new WebService(GIP_HTTP)
    def meta = [:];
    meta.samples = args.collect { sampleId -> 

        println "Query Assets for $sampleId ..."

        def giphubSampleId = sampleId + 'RNA'
        def result = ws.path("xip/sample/${giphubSampleId}/fastq").get(runid:'all')

        if(result == null)
            throw new bpipe.PipelineError("Unablet o query fastq for sample id $giphubSampleId")

        // export API is the only one that will resolve lims id without family id
        def assetResult = ws.path("dashboard/export/resolve")
                            .post(bam:true,lims_ids:sampleId,autoselect:true)
                            .assets.groupBy { it.filetype }


        println "Response: " + assetResult

        def limsId_familyId = assetResult.bam[0].assets[0].metadata.lims_id

        println "LimsId-FamilyId for $sampleId is $limsId_familyId"

        assert limsId_familyId.startsWith(sampleId + '-')

        def familyId = limsId_familyId.tokenize('-')[1]

        def sex = ws.path("xip/sample/$limsId_familyId/sex").get()[limsId_familyId].toLowerCase()

        def fastqs = result[giphubSampleId];
        def bam = assetResult.bam[0].assets[0].fullpath
        def rare_vcf = ws.path("xip/sample/$sampleId-$familyId/asset?type=vcf").get().path
        def full_vcf = ws.path("xip/sample/$sampleId-$familyId/asset?type=norm.vcf.gz").get().path
        def pedigree = ws.path("xip/sample/$sampleId-$familyId/pedigree").get()
        if(pedigree)
            pedigree = pedigree[limsId_familyId]

        def sampleRes = [
            identifier : sampleId,
            familyId : familyId,
            sex : sex,
            fastqs : fastqs,
            wgs_bam : bam,
            rare_vcf : rare_vcf,
            full_vcf : full_vcf,
            parents  : pedigree?.tokenize('=')?.getAt(1)?.tokenize(',')?.collect { it.trim().tokenize('-')[0] }
        ]
        return sampleRes;
    }

    if(!CONTROL_SET) {
        throw new bpipe.PipelineError("CONTROL_SET variable was not set. Please configure this variable in config.groovy to point to the directory containing control sample sheets")
    }

    meta.control_groups = [
        adult_male: "$CONTROL_SET/sample-sheet-adult_male.tsv".toString(),
        adult_female: "$CONTROL_SET/sample-sheet-adult_female.tsv".toString(),
        child_male: "$CONTROL_SET/sample-sheet-child_male.tsv".toString(),
        child_female: "$CONTROL_SET/sample-sheet-child_female.tsv".toString(),
    ]

    new File('samples.yaml').withWriter { writer -> new Yaml().dump(meta, writer) };
    return readMetaFromFile('samples.yaml')
}

