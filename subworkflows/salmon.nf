import java.nio.file.Files

include { SALMON_INDEX } from '../modules/salmon/index'
include { SALMON_QUANT } from '../modules/salmon/quant'
include { TAR as TAR_SALMON_INDEX } from '../modules/utils/tar'
include { TAR as TAR_SALMON_QUANT } from '../modules/utils/tar'

workflow SALMON_RNASEQ {
    take:
        samples          // Channel with tuples of the form: [val(sample_meta), [read1, read2]] or [val(sample_meta), [read1]]
        transcript_fasta // NULL or Channel: /path/to/transcriptome/fasta
        genome_fasta     // NULL or Channel: /path/to/genome/fasta
        index_path       // NULL or string pointing to: /path/pre-build/salmon/index

    main:
        // Find the index or build it from the FASTA files
        if (index_path == "" || index_path == null || file(index_path).exists() == false) {
            
            // build
            index = SALMON_INDEX(transcript_fasta, genome_fasta)
            salmon_index = index.index
            TAR_SALMON_INDEX(salmon_index)

            // Quantify with Salmon
            quant = SALMON_QUANT(samples, salmon_index)

        } else {

            // Find index
            salmon_index = Channel.fromPath(index_path)

            // Use first() to keep repeating the output. This workaround is needed in case the user passed the 
            // index file manually. If absent, nextflow would then only pass it to the first item. By using first()
            // we can re-use the item
            quant = SALMON_QUANT(samples, salmon_index.first())
        }

        // compress output (directory paths only)
        TAR_SALMON_QUANT(quant.results.map{ it[1] })

    emit:
        index   = salmon_index
        results = quant.results
        sams    = quant.sams
}
