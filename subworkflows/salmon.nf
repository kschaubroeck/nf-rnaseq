import java.nio.file.Files

include { SALMON_INDEX } from '../modules/salmon/index'
include { SALMON_QUANT } from '../modules/salmon/quant'
include { GUNZIP as GUNZIP_INDEX } from '../modules/utils/gunzip'

workflow SALMON_RNASEQ {
    take:
        samples          // Channel with tuples of the form: [val(sample_meta), [read1, read2]] or [val(sample_meta), [read1]]
        transcript_fasta // NULL or Channel: /path/to/transcriptome/fasta
        genome_fasta     // NULL or Channel: /path/to/genome/fasta
        index_path       // NULL or string pointing to: /path/pre-build/salmon/index

    main:
        // Find the index or build it from the FASTA files
        if (index_path == "" || index_path == null || file(index_path).exists() == false) {
            salmon_index = SALMON_INDEX(transcript_fasta, genome_fasta)
        } else {
            salmon_index = Channel.fromPath(index_path)

            // If the index is gziped, then ungzip
            if (Files.probeContentType(salmon_index) == "application/gzip") { index = GUNZIP_INDEX(index) }
        }

        // Quantify with Salmon
        quant = SALMON_QUANT(samples, salmon_index.index)

    emit:
        index   = salmon_index.index
        results = quant.results
        sams    = quant.sams
}
