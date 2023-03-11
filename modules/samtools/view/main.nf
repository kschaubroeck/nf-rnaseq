process SAMTOOLS_VIEW {
    tag   "$sample"
    label "process_medium_resources"
    label "tool_samtools"

    input:
        val  sample
        path sam_file

    output:
        tuple val(sample), path(bam), emit: bams

    script:
        args = task.ext.args ?: ""
        bam  = task.ext.bamFile ?: sam_file.getSimpleName() + ".bam"
        """
        samtools view $args -@$task.cpus --with-header --bam --output $bam $sam_file
        """
    stub:
        bam  = task.ext.bamFile ?: sam_file.getSimpleName() + ".bam"
        """
        touch $bam
        """
}
