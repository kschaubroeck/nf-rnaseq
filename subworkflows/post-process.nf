import java.nio.file.Files

include { TERMINUS_GROUP } from '../modules/terminus/group'
include { TERMINUS_COLLAPSE } from '../modules/terminus/collapse'
include { TXIMETA_FISHPOND } from '../modules/R/tximeta-fishpond'

workflow POST_PROCESS_SALMON_RNASEQ {
    take:
        salmon  // results cahnnel from SALMON_QUANT

    main:
        directories = salmon.map { it[1] }.collect()
        samples = salmon.map { it[0].id }.collect()

        // Begin tximeta and fishpond
        TXIMETA_FISHPOND(samples, directories)

        // Terminus for grouping together problematic transcripts
        groups = TERMINUS_GROUP(salmon)
        TERMINUS_COLLAPSE(groups.salmon.collect(), groups.terminus.collect())
}

