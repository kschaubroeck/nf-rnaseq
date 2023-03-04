process GUNZIP {
    tag "$input"
    label "ubuntu"

    input:
        path input

    output:
        path outname

    script:
        args    = task.ext.args ?: ""

        // Predict output name by removing the .gz extnesion
        outname = input.lastIndexOf('.') != -1 ? input.substring(0, input.lastIndexOf(".")) : input
        """
        gunzip $args $input
        """

    stub:
        outname = input.lastIndexOf('.') != -1 ? input.substring(0, input.lastIndexOf(".")) : input
        """
        touch $outname
        """
}
