process SALMON_QUANT {
    tag "$sample.id"
    label "process_medium_resources"
    label "tool_salmon"

    input:
        tuple val(sample), path(reads)
        path  index

    output:
        tuple val(sample), path("$sample.id"), emit: results
        stdout                                 emit: sams

    script:
        args   = task.ext.args ?: ""
        single = (sample.single == true ? "--fldMean $sample.mean --fldSD $sample.sd" : "")
        fastqs = (sample.single == true ? "--unmatedReads ${reads[0]}" : "--mates1 ${reads[0]} --mates2 ${reads[1]}")
        """
        salmon quant $args --dumpEqWeights --threads $task.cpus --index $index $fastqs $single --output $sample.id
        """

    stub:
        """
        mkdir $sample.id
        touch $sample.id/quant.sf
        echo 'CONTENT FOR SAM FILE'
        """
}
