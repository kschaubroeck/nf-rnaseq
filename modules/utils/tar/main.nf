process TAR {
    tag "$input"
    label "ubuntu"

    input:
        path input

    output:
        path outname

    script:
        args    = task.ext.args ?: ""
        outname = input.getBaseName() + ".tar"
        """
        tar -czvf $args $outname $input
        """

    stub:
        outname = input.getBaseName() + ".tar"
        """
        touch $outname
        """
}
