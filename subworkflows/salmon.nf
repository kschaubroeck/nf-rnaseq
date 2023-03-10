import java.nio.file.Files

include { SALMON_INDEX } from '../modules/salmon/index'
include { SALMON_QUANT } from '../modules/salmon/quant'

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

    emit:
        index   = salmon_index
        results = quant.results
        sams    = quant.sams
}
