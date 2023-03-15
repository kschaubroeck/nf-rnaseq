process TXIMETA_FISHPOND {
    tag "$samples"
    label "process_high_resources"
    label "tool_tximeta_fishpond"

    input:
        val  samples
        path paths

    output:
        path "transcript/*",                      emit: transcripts
        path "gene/*",                            emit: gene
        path "dtu/*",                             emit: dtu
        path "quantification-summary.csv.gz",     emit: quant
        path "library-format.csv.gz",             emit: library_format
        path "unique-mappings.csv.gz",            emit: unique_mappings
        path "ambiguous-mappings.csv.gz",         emit: ambiguous_mappings
        path "tximeta.log",                       emit: log
        path "tximeta.json",                      emit: tximeta_json
        path "annotation/*",                      emit: annotation_data
        path "sequence/*",                        emit: sequence
        path "*.sqlite",                          emit: db

    script:
        samples = samples.join(' ')
        paths   = paths.join(' ')
        args   = task.ext.args ?: ""
        """
        mkdir -p transcript
        mkdir -p dtu
        mkdir -p gene
        mkdir -p annotation
        mkdir -p annotation/tsv
        mkdir -p annotation/rds
        mkdir -p sequence
        
        tximeta-fishpond.r --samples $samples --directory $paths $args &> tximeta.log
        """

    stub:
        samples = samples.join(' ')
        paths   = paths.join(' ')
        """
        mkdir -p transcript
        mkdir -p dtu
        mkdir -p gene
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

        touch transcript/effective-lengths.csv.gz
        touch gene/effective-lengths.csv.gz
        touch dtu/effective-lengths.csv.gz

        touch transcript/infreps-mean.csv.gz
        touch gene/infreps-mean.csv.gz
        touch dtu/infreps-mean.csv.gz

        touch transcript/infreps-relative-variance.csv.gz
        touch gene/infreps-relative-variance.csv.gzeffective-lengths.csv.gz
        touch dtu/infreps-relative-variance.csv.gzeffective-lengths.csv.gz

        touch transcript/raw-counts.csv
        touch gene/raw-counts.csv

        touch transcript/scaled-counts.csv
        touch gene/scaled-counts.csv
        touch dtu/scaled-counts.csv

        touch transcript/infreps.rds
        touch gene/infreps.rds
        touch dtu/infreps.rds

        touch quantification-summary.csv.gz
        touch library-format.csv.gz
        touch unique-mappings.csv.gz
        touch ambiguous-mappings.csv.gz

        stub-tximeta-fishpond.r --samples $samples --directory $paths &> tximeta.log
        """
}
