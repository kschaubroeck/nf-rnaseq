process TERMINUS_GROUP {
    tag "$sample.id"
    label "process_medium_resources"
    label "tool_terminus"

    input:
        tuple val(sample), path(quant)

    output:
        path salmon_output,     type: "dir", emit: salmon
        path outputDir,         type: "dir", emit: terminus

    script:
        args = task.ext.args ?: ""

        base = quant.getBaseName()
        outputDir = "terminus_output/$base"

        salmon_output = quant
        """
        terminus group $args -d $salmon_output -o terminus_output
        """

    stub:
        salmon_output = quant
        base = quant.getBaseName()
        outputDir = "terminus_output/$base"
        """
        mkdir terminus_output
        mkdir $outputDir
        touch $outputDir/groups.txt
        """
}
