process TXIMETA_FISHPOND {
    tag "$samples"
    label "process_medium_resources"
    label "tool_tximeta_fishpond"

    input:
        val  samples
        path paths

    output:
        path "transcript/*",                      emit: transcripts
        path "gene/*",                            emit: gene
        path "quantification-summary.csv.gz",     emit: quant
        path "library-format.csv.gz",             emit: library_format
        path "unique-mappings-summary.csv.gz",    emit: unique_mappings
        path "ambiguous-mappings-summary.csv.gz", emit: ambiguous_mappings
        path "tximeta.log",                       emit: log
        path "tximeta.json",                      emit: tximeta_json
        path "annotation/*",                      emit: annotation_data
        path "sequence/*",                        emit: sequence
        path "*.sqlite",                          emit: db

    script:
        samples = samples.join(' ')
        paths   = paths.join(' ')
        """
        tximeta-fishpond.r --samples $samples --directory $paths &> tximeta.log
        """

    stub:
        samples = samples.join(' ')
        paths   = paths.join(' ')
        """
        mkdir -p transcript
        mkdir -p transcript/infreps
        mkdir -p gene
        mkdir -p gene/infreps
        mkdir -p annotation
        mkdir -p annotation/tsv
        mkdir -p annotation/rds
        mkdir -p sequence

        touch annotation/metadata.json
        touch annotation/source.json
        touch tximeta.json

        touch annotation/tsv/transcripts.tsv.gz
        touch annotation/rds/transcripts.rds
        touch annotation/rds/genes.rds
        touch annotation/tsv/genes.tsv.gz

        touch sequence/cdna.rds
        touch SOURCE.v000.GrCh38.sqlite

        touch transcript/counts-length-scaled.csv
        touch transcript/abundance-tpm.csv
        touch transcript/effective-lengths.csv
        touch transcript/infreps/raw-variance.csv
        touch transcript/infreps/counts-length-scaled.rds

        touch transcript/effective-lengths.csv.gz
        touch gene/effective-lengths.csv.gz

        touch transcript/infreps-mean.csv.gz
        touch gene/infreps-mean.csv.gz

        touch transcript/infreps-relative-variance.csv.gz
        touch gene/infreps-relative-variance.csv.gzeffective-lengths.csv.gz

        touch transcript/scaled-counts.csv
        touch gene/scaled-counts.cs

        touch transcript/infreps.rds
        touch gene/infreps.rds

        touch quantification-summary.csv.gz
        touch library-format.csv.gz
        touch unique-mappings-summary.csv.gz
        touch ambiguous-mappings-summary.csv.gz

        stub-tximeta-fishpond.r --samples $samples --directory $paths &> tximeta.log
        """
}
