process GZIP {
    tag "$input"
    label "ubuntu"

    input:
        path input

    output:
        path outname

    script:
        args    = task.ext.args ?: ""
        outname = input + ".gz"
        """
        gzip $args $input
        """

    stub:
        outname = input + ".gz"
        """
        touch $outname
        """
}
