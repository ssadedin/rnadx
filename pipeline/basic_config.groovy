load "config_rnadx_home.groovy"

STAR_GENOME_DIR="$RNADX_BASE/reference/hg38/star/star_2.7.3a/" 
GENOME_FASTA="$RNADX_BASE/reference/hg38/genome/GRCh38.primary_assembly.genome.fa"
GENOME_ANNOTATION="$RNADX_BASE/reference/hg38/annotation/gencode.v37.annotation.gtf"

//
//trimmomatic parameters
//

adapters = "$RNADX_BASE/reference/trimmomatic/adapters/NexteraPE-PE.fa"
LEADING = 25
TRAILING = 25
MINLEN = 50
ILLUMINACLIP = "2:30:10"

//
// Multiqc parameters
//

MULTIQCDIR = "."

//
// STAR parameters
//

// length of genomic sequence around annotated junctioin; max(ReadLength)-1; default: 100
SJBOHANG = 99

//
// DESeq2 parameters
//

DESIMPLE_HOME = "$RNADX_BASE/deseq2-simple"

//
// Global parameters
//

//number of treads per sample for multithreaded tools
NTHREADS = 4

