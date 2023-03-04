process SAMTOOLS_VIEW {
    tag   "$sample"
    label "process_medium_resources"
    label "tool_samtools"

    input:
        val sample
        stdin

    output:
        tuple val(sample), path(bam), emit: bams

    script:
        args = task.ext.args ?: ""
        bam  = task.ext.bamFile ?: sample + ".bam"
        """
        cat - | samtools view $args -@$task.cpus --bam --output $bam --with-header
        """

    stub:
        bam  = task.ext.bamFile ?: sample + ".bam"
        """
        cat - > $bam
        """
}
