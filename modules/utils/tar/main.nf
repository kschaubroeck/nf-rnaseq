process TAR {
    tag "$input"
    label "ubuntu"

    input:
        path input

    output:
        path outname, emit: compressed

    script:
        args    = task.ext.args ?: ""
        outname = task.ext.fileName ?: (input.getName() + ".tar")
        """
        tar -czvf $args $outname $input
        """

    stub:
        outname = task.ext.fileName ?: (input.getName() + ".tar")
        """
        touch $outname
        """
}
