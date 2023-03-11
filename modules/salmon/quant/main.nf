process SALMON_QUANT {
    tag "$sample.id"
    label "process_medium_resources"
    label "tool_salmon"

    input:
        tuple val(sample), path(reads)
        path  index

    output:
        tuple val(sample), path("$sample.id"), emit: results
        path  sam_file,                        optional: true, emit: sams

    script:
        args   = task.ext.args ?: ""
        single = (sample.single == true ? "--fldMean $sample.mean --fldSD $sample.sd" : "")
        fastqs = (sample.single == true ? "--unmatedReads ${reads[0]}" : "--mates1 ${reads[0]} --mates2 ${reads[1]}")
        sam_file = sample.id + ".sam.gz"

        if (task.ext.outSam == true) {
            """
            salmon quant --writeMappings $args --dumpEqWeights --threads $task.cpus --index $index $fastqs $single --output $sample.id | gzip -9 > $sam_file
            """           
        } else {
            """
            salmon quant $args --dumpEqWeights --threads $task.cpus --index $index $fastqs $single --output $sample.id 
            """
        }
    stub:
        sam_file = sample.id + ".sam.gz"
        if (task.ext.outSam == true) {
            """
            mkdir $sample.id
            touch $sample.id/quant.sf
            touch $sam_file
            """       
        } else {
            """
            mkdir $sample.id
            touch $sample.id/quant.sf
            """  
        }
}
