process TAR {
    tag "$input"
    label "process_medium_resources"
    label "ubuntu"

    input:
        path input

    output:
        path outname, emit: compressed

    script:
        args    = task.ext.args ?: ""
        outname = task.ext.fileName ?: (input.getName() + ".tar.gz")
        """
        tar $args -czhf $outname $input
        """

    stub:
        outname = task.ext.fileName ?: (input.getName() + ".tar.gz")
        """
        touch $outname
        """
}
