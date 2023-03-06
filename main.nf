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
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SALMON_ALIGNMENT } from './modules/samtools/view'

/* Subworkflow includes
---------------------------------------------------------------------------- */
include { RNASEQ_PRE_PROCESS } from './subworkflows/pre-process'
include { SALMON_RNASEQ } from './subworkflows/salmon'
include { POST_PROCESS_SALMON_RNASEQ } from './subworkflows/post-process'

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
        trackDiscardedReads: "Write a copy of all reads thrown out due to filtering and track which reads became orphaned.",
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
    // Load the samples from CSV file
    samples = RNASEQ_PRE_PROCESS(params.csv)

    // Align and quantify with Salmon
    salmon  = SALMON_RNASEQ(samples.trimmed, params.transcriptome, params.genome, params.index)

    // Salmon alignments to BAM
    if (params.bam) {
        SAMTOOLS_VIEW_SALMON_ALIGNMENT(salmon.results.map{ it[0].id }, salmon.sams)
    }

    // Post-process quants
    POST_PROCESS_SALMON_RNASEQ(salmon.results)
}

