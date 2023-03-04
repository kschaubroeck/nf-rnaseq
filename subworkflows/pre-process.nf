import java.nio.file.Files

include { import_samples } from '../modules/samples/import'
include { FASTP } from '../modules/fastp'

workflow RNASEQ_PRE_PROCESS {
    take:
        samples_csv // Channel: csv file of samples

    main:
        // Make sure the sample csv file exists
        if (samples_csv == "" || file(samples_csv).exists() == false) {
            error "Sample CSV file located at $samples_csv does not exist."
        }

        // Import load samples and pre-process using fastp
        fastp = FASTP(import_samples(samples_csv))
    emit:
        trimmed = fastp.trimmed
        json    = fastp.json
        html    = fastp.html
}
