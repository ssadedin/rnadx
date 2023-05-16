import org.yaml.snakeyaml.*

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

def image_list = [
        "dockanomics/essential-r:4.0",
        "dockanomics/fraser:1.2",
        "dockanomics/star:2.7",
        "dockanomics/scindo:latest",
        "maplesond/portcullis:latest"
    ]

set_image = {
    branch.image = branch.name
}

fetch_singularity_image = {
    output.dir = "singularity"
    def mangled = image.replaceAll("[-/:.]", "_")
    produce(mangled + ".sif") {
        exec """
            singularity pull --disable-cache $output docker://$image
        ""","singpull"
    }
}

fetch_mendeliome_panel = {
    output.dir = "reference/panels"
    def mangled = image.replaceAll("[-/:.]", "_")
    produce("Mendeliome.tsv") {
        exec """
            curl -s -o $output https://panelapp.agha.umccr.org/panels/137/download/01234/
        ""","singpull"
    }
}

install_fraser = {
    exec """
        export R_LIBS_USER=$RNADX_BASE/tools/r_libs

        mkdir -p $R_LIBS_USER

        Rscript --vanilla -e 'BiocManager::install("FRASER")'
    ""","singpull"
}

run() {
    fetch_mendeliome_panel +
    image_list * [set_image + fetch_singularity_image] +
    star_genome_gen +
    install_fraser
}
