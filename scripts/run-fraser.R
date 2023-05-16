library(data.table)
library(ggplot2)
library(FRASER)
library(GenomicFeatures)
library(R.utils)

batch.name <- Sys.getenv("RNADX_BATCH");
gencode.annotation.gtf <- Sys.getenv("GENCODE_ANNOTATION_GTF");
control.sheet  <- Sys.getenv("FRASER_CONTROLS");
test.sheet  <- Sys.getenv("FRASER_TEST");
output.file  <- Sys.getenv("FRASER_OUTPUT_FILE");
working.dir <- Sys.getenv("FRASER_WORKING_DIR");

#batch.base <- paste0("/batches/", batch.name);
batch.base <- Sys.getenv("BATCH_BASE");

msg <- function(m,...) {
    print(sprintf(m,...))
}

setwd(batch.base);

options(width=1000)

load_sheet = function(sheet) {
    msg("Loading sheet %s", sheet)
    batch_dir = dirname(dirname(sheet))
    msg("Batch dir is %s", batch_dir)
    sampleTable <- fread(sheet, header=T)
    sampleTable[, sampleID := sample]
    #SS: getAbsolutePath is returning incorrect path combo of control file + working dir
    #sampleTable[, bamFile := getAbsolutePath(bam, workingDirectory=batch_dir)];
    #sampleTable[, bamFile := paste0(batch_dir, "/", bam)];
    sampleTable$pairedEnd <- TRUE
    return(sampleTable)
}

controlTable = load_sheet(control.sheet)
testSampleTable = load_sheet(test.sheet)

msg("Loaded sample sheets") 

sampleTable = rbind(testSampleTable,controlTable)
sampleTable <- sampleTable[!grepl('NTC', sampleTable$bam),]
sampleTable$group = 1
sampleTable$bamFile = sampleTable$bam

msg("Constructed sample table")

print(sampleTable)

full.working.dir=paste0(batch.base, "/fraser/", working.dir)

msg("Creating FraserDataSet with working dir %s", full.working.dir)

settings <- FraserDataSet(colData=sampleTable, workingDir=full.working.dir);

msg("Constructed fraser data set")

# deliberately setting this higher than 
# the core count we are using because empirically
# it spends a lot of time under-utilising the cores
# as some samples process faster than others
# from inspecting the code many stages are broken down
# by chromosome or cores passed through to subread so this
# can let us get more usage from the pool of cores assigned.
# need to monitor and make sure we are not over-subscribing and
# losing out but based on current behaviour it looks unlikely
register(MulticoreParam(workers=24))

msg("Starting to count RNA data")

fds <- countRNAData(settings);

fds <- calculatePSIValues(fds)
fds <- filterExpressionAndVariability(fds, minDeltaPsi=0.1, filter=TRUE)
fds <- fit(fds, q=3, type="psi5", implementation="PCA", iterations=2)
fds <- calculatePvalues(fds, type="psi5")
fds <- calculatePadjValues(fds, type="psi5", method="BY")
fds <- calculateZscore(fds, type="psi5")

res <- results(fds, psiType='psi5', zScoreCutoff=NA, padjCutoff=NA, deltaPsiCutoff=NA);

# Read in Gencode, and extract just the gene entries.
#
genc <- rtracklayer::import(gencode.annotation.gtf);
genc.genes <- genc[genc$type == "gene"]

# Find genes which overlap the results,
# and find which results had an overlapping gene
#
annot.hits <- findOverlaps(res, genc.genes, type="any", select="arbitrary")
annot.rows <- which(!is.na(annot.hits))

# Make a table of gene rows corresponding to the
# results we want to annotate.
#
genc.rows <- genc.genes[annot.hits[annot.rows]]

# Copy across the strand information
#
annot.strand <- strand(res);
annot.strand[annot.rows] <- strand(genc.rows)
strand(res) <- annot.strand;

# Copy across the gene name.
#
res$gene_name <- NA;
res[annot.rows]$gene_name <- genc.rows$gene_name;

# Now make a table of known junctions
#
genc.exons <- as.data.table(genc[genc$type=="exon"])
known.junctions <- genc.exons[order(gene_id, start, end),
                              .(seqnames = head(seqnames, -1),
                                start = head(end, -1) + 1,
                                end = tail(start, -1) - 1,
                                strand = head(strand, -1)),
                              by=.(transcript_id)];
known.junctions <- unique(known.junctions[, .(seqnames, start, end, strand, known = TRUE)])

res.tbl <- as.data.table(res);
res.tbl <- merge(res.tbl, known.junctions, by=c("seqnames", "start", "end", "strand"), all.x=TRUE);
res.tbl[is.na(known), known := FALSE];
res.tbl.cols <- colnames(res.tbl);
res.tbl.cols <- res.tbl.cols[res.tbl.cols != "hgncSymbol" & res.tbl.cols != "type"];
res.tbl <- res.tbl[, ..res.tbl.cols]
res.tbl <- res.tbl[order(padjust, pValue, -counts)]
fwrite(res.tbl, paste0(output.file, ".nonsig"), sep="\t")

significant.results = res.tbl[abs(res.tbl$zScore)>1.5 & res.tbl$pValue<0.05 & (res.tbl$sampleID %in% testSampleTable$sample) & (abs(res.tbl$deltaPsi) / res.tbl$psiValue)>0.1 ,]

fwrite(significant.results, output.file, sep="\t")

