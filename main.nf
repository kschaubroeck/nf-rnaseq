#! usr/bin/env nextflow

/* ANSI color Codes
---------------------------------------------------------------------------- */
// Used for determining the color of text output to the console
ANSI_RESET  = "\u001B[0m"
ANSI_GREEN  = "\u001B[32m"
ANSI_RED    = "\u001B[31m"
ANSI_YELLOW = "\u001B[33m"
ANSI_BLUE   = "\u001B[34m"
ANSI_GRAY   = "\033[0;90m"
ANSI_GREY   = ANSI_GRAY

/* Module includes
---------------------------------------------------------------------------- */
include { import_samples } from './modules/samples/import'
include { FASTP } from './modules/fastp'
include { SALMON_INDEX } from './modules/salmon/index'
include { SALMON_QUANT } from './modules/salmon/quant'
include { TAR as TAR_SALMON_INDEX } from './modules/utils/tar'
include { TAR as TAR_SALMON_QUANT } from './modules/utils/tar'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SALMON_ALIGNMENT } from './modules/samtools/view'
include { TERMINUS_GROUP } from './modules/terminus/group'
include { TERMINUS_COLLAPSE } from './modules/terminus/collapse'
include { TXIMETA_FISHPOND } from './modules/R/tximeta-fishpond'
include { TAR as TAR_TERMINUS } from './modules/utils/tar'

/* Help and start-up
---------------------------------------------------------------------------- */
log.info"""
██████╗░███╗░░██╗░█████╗░░░░░░░░██████╗███████╗░██████╗░
██╔══██╗████╗░██║██╔══██╗░░░░░░██╔════╝██╔════╝██╔═══██╗
██████╔╝██╔██╗██║███████║█████╗╚█████╗░█████╗░░██║██╗██║
██╔══██╗██║╚████║██╔══██║╚════╝░╚═══██╗██╔══╝░░╚██████╔╝
██║░░██║██║░╚███║██║░░██║░░░░░░██████╔╝███████╗░╚═██╔═╝░
╚═╝░░╚═╝╚═╝░░╚══╝╚═╝░░╚═╝░░░░░░╚═════╝░╚══════╝░░░╚═╝░░░

${ANSI_GREEN}${workflow.manifest.name}${ANSI_RESET} version ${ANSI_YELLOW}${workflow.manifest.version}${ANSI_RESET} by ${ANSI_BLUE}${workflow.manifest.author}${ANSI_RESET}
${ANSI_GREEN}Script:${ANSI_RESET} ${ANSI_YELLOW}${workflow.scriptName}${ANSI_RESET}
"""
/* Help code
---------------------------------------------------------------------------- */
help_messages = [
    "Global": [
        help: "Display this help message",
        csv: "A CSV file containing a list of samples and their FASTQ files."
    ],
    "Fastp": [
        skipFastqPrep: "Skip pre-processing FASTQ files using fastp.",
        saveTrimmedFastq: "Save a copy of the processed FASTQ files.",
        trimmedFastqCompression: "Level of compression to use when gzipping FASTQ files",
        baseQuality: "Minimum phred score as base must achieve to 'pass' this filte.",
        requiredQualityBases: "Percent of bases that must pass the quality filter to keep the read",
        unknownBaseLimit: "Maximum number of N bases a read can have or else it will be thrown out",
        disableBaseQualityFilter: "Disable quality inspection of bases",
        defineOverlapLength: "Amount of nucleotides that defines an overlap",
        defineOverlapDifferenceNumber: "Maximum number of bases that can be different in an overlap",
        defineOverlapPercentDiffernce: "Maximum percent of bases that can be different in an overlap",
        enableBaseCorrectionByOverlap: "Enable replacement of low-quality bases if their is an overlapping complement of highn quality.",
        autoDetectAdapters: "Enable auto-detection of adapters based on read overlap.",
        adapterFasta: "File containing the FASTA sequences of adapters.",
        adapterSequence: "Sequence of the adapter on the first (R1) read.",
        adapterSequence2: "Sequence of the adapter on the second (R2) read.",
        disableAdapterTrimming: "Disable adapter trimming",
        removeDuplicates: "Remove duplicate reads",
        disableDuplicationEval: "Do not test of duplicate reads.",
        duplicateSearchAccuracy: "Accuracy level on a 1-6 scale. Higher levels are more accurate but use more memory.",
        enableLengthFilter: "Enable filtering of reads based on length",
        minReadLength: "Minimum size of a read",
        maxReadLength: "Maximum size of a read",
        disableRepresentedSeqAnalysis: "Disable testing for over-representation of sequences.",
        representedSeqSamplingRate: "Frequency to take reads and sample their sequences for representation testing"
    ],
    "Index": [
        genome: "Path or URL to genome FASTA file",
        transcriptome: "Path or URL to transcript cDNA sequences",
        kmer: "Length of k-mer sequences when building salmon-index",
        gencode: "States the cDNA FASTA files came from gencode and therefore need extra processing.",
        index: "Path to a pre-built Salmon index"
    ],
    "Alignment": [
        single: "Flag to set reads as single-end",
        mean: "Mean size of a fragment in the library",
        sd: "Standard deviation of fragment sizes",
        libtype: "Library strandedness and read orientation",
        allowDovetail: "Allow dovetailed reads. Dovetailed reads are reads that extend past eachother.",
        disableOrphanRecovery: "Disable recovery of orphaned reads.",
        salmonClip: "Lightly clip the ends of reads so mismatches do not count as penalties",
        salmonClipOverhangs: "Clip overhangs that may belong to several transcripts to avoid penalties",
        mismatchSkipNumber: "When extending k-mers and a mismatch is encountered, skip a few bases to continue searching.",
        seqSimilarityThreshold: "Maximum number of places a MEM can appear in the genome",
        readOccuranceThreshold: "Maximum number of places a read can map to",
        noChainingLimit: "Disable the limited number of iterations when combining MEM chains.",
        minQualityScore: "Quality score a read or read-pair must achieve to be counted",
        bam: "Write the mappings out as BAM files.",
        writeUnmapped: "Write a text files that contains names of reads that did not map anywhere.",
        unmappedWithQualityStrings: "Write the quality strings when outputting the BAM mappings.",
        writeOrphans: "Write the transcripts who contain orphaned links."
    ],
    "Quant": [
        disableGCBiasCorrection: "Disable correcting for GC biases.",
        disableHexamerBiasCorrection: "Disable biases correcting for random hexamer (position) priming",
        disablePositionBiasCorrection: "Disable bias correction for fragment locations",
        incompatProbability: "Probability that a fragmet mapping is correct even if the orientation/strandedness is wrong",
        factorization: "Number of bins to factor transcript equivalnce classes into",
        gibbs: "Number of gibbs samples of the posterior to take",
        disableGammaDraw: "Disable use of the gamma distribution when estimating counts from the posterior.",
        thinningFactor: "Number of iterations to skip before recording a gibbs sampling.",
        bootstraps: "Number of bootstraps to take"
    ],
    "Terminus": [
        groupMinPosteriorSpread: "Spread of posterior a transcript must achieve to be considered for grouping",
        terminusWeighttolerance: "Difference in conditional probability when defining it two groups are similar",
        terminusSeed: "Random seed",
        minSamplesForCollapse: "Percent of samples a group must occur in to be considered for collapsing"
    ],
    "Fishpond-Tximeta": [
        fishpondMinCount: "Minimum count a transcript/gene must get to be included in a size factor calculation",
        fishpondMinN: "Minimum number of samples that must hit the minimum count for a gene/transcript to be used in size calculations",
        fishpondMinCountGene: "fishpondMinCount, but for genes",
        fishpondMinNGene: "fishpondMinNGene, but for genes"
    ]
]

if (params.help) {
    include { show_help } from './modules/utils/help'
    show_help(help_messages)
    exit 0
}

/* Log the selected parameters
---------------------------------------------------------------------------- */
include { log_params } from './modules/utils/log'
log_params(help_messages, params)

/* Main Workflow
---------------------------------------------------------------------------- */
workflow {
    // Make sure the sample csv file exists
    if (params.csv == "" || file(params.csv).exists() == false) {
        error "Sample CSV file located at $params.csv does not exist."
    }

    // Import load samples and pre-process using fastp
    samples = FASTP(import_samples(params.csv))

    // Find the index or build it from the FASTA files
    if (params.index == "" || params.index == null || file(params.index).exists() == false) {
        // Build index
        salmon_index = SALMON_INDEX(params.transcriptome, params.genome)
        TAR_SALMON_INDEX(salmon_index.index)

        // Quantify
        quant = SALMON_QUANT(samples.trimmed, salmon_index.index)
    } else {
        // Pre-supplied index
        salmon_index = Channel.fromPath(params.index)

        // Use first() to keep repeating the output. This workaround is needed in case the user passed the 
        // index file manually. If absent, nextflow would then only pass it to the first item. By using first()
       // we can re-use the item
        quant = SALMON_QUANT(samples.trimmed, salmon_index.first())
    }

    // compress output (directory paths only)
    TAR_SALMON_QUANT(quant.results.map{ it[1] })

    // Salmon alignments to BAM
    if (params.bam) {
        SAMTOOLS_VIEW_SALMON_ALIGNMENT(quant.results.map{ it[0].id }, quant.sams)
    }

    // Prepare for post-processing
    directories = quant.results.map { it[1] }.collect()
    samples = quant.results.map { it[0].id }.collect()

    // Begin tximeta and fishpond
    TXIMETA_FISHPOND(samples, directories)

    // Terminus for grouping together problematic transcripts
    groups = TERMINUS_GROUP(quant.results)
    terminus = TERMINUS_COLLAPSE(groups.salmon.collect(), groups.terminus.collect())

    // compress terminus output to save space
    // We need to get all child items, which is why we flatten the channel
    TAR_TERMINUS(terminus.results.flatten())
}

