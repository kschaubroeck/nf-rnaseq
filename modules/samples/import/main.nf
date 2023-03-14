
def import_samples(samples_metadata, seperator = ",") {
    // Fetch samples from file
    samples = Channel.fromPath(samples_metadata, checkIfExists: true)
                .splitCsv(sep: seperator, header: ["sampleId", "read1", "read2", "mean", "sd"])
                .take( params.dev ? params.devInputs : -1 )
    .map { row ->
        // Let's get the error checking out of the way and make sure there is an id
        if (row.sampleId == null || row.sampleId == "") {
            error "Error importing samples: At least one sample is missing an ID code."
        }

        // Remove spaces in names
        if (row.sampleId.indexOf(" ") != -1) {
            log.info "WARNING: Spaces were detected for sample '$row.sampleId' and will be replaced."
            row.sampleId = row.sampleId.replace(" ", "_")
        }

        // Now let's get checking the first read
        if (row.read1 == null || row.read1 == "") {
            error "Error importing samples: Read1 is undefined for $row.sampleI.d"
        }

        // Get a file and check fome properties
        read1 = file(row.read1)

        if (read1.exists() == false) {
            error "Error importing samples: Read1 does not exist for $row.sampleId."
        }

        if (read1.getName().indexOf(".fastq") == -1 && read1.getName().indexOf(".fa") == -1) {
            error "Error importing samples: Read1 from $row.sampleId does not contain the fastq extension."
        }

        // Now do some checks on read2
        read2 = ""
        if (row.read2 != null && row.read2 != "") {
            read2 = file(row.read2)

            if (read2.exists() == false) {
                error "Error importing samples: Read2 does not exist for $row.sampleId."
            }

            if (read2.getName().indexOf(".fastq") == -1 && read1.getName().indexOf(".fa") == -1) {
                error "Error importing samples: Read2 from $row.sampleId does not contain the fastq extension."
            }
        }

        // Initialize variables and populate the data with the sample's ID
        def meta = [:]
        meta.id = row.sampleId
        meta.single = (row.read1 != "" && row.read1 != null && (row.read2 == "" || row.read2 == null)) ? true : false
        meta.mean     = row.mean ?: null
        meta.sd       = row.sd   ?: null

        // // Return in the appropriate forms
        if (meta.single) {
            return [meta, [read1]]
        }

        return [meta, [read1, read2]]
    }
}
