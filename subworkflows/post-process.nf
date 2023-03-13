import java.nio.file.Files

include { TERMINUS_GROUP } from '../modules/terminus/group'
include { TERMINUS_COLLAPSE } from '../modules/terminus/collapse'
include { TXIMETA_FISHPOND } from '../modules/R/tximeta-fishpond'
include { TAR as TAR_TERMINUS } from '../modules/utils/tar'

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
        terminus = TERMINUS_COLLAPSE(groups.salmon.collect(), groups.terminus.collect())

        // compress terminus output to save space
        // We need to get all child items, which is why we flatten the channel
        TAR_TERMINUS(terminus.results.flatten())
}

