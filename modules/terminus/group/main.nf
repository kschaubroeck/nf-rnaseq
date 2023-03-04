process TERMINUS_GROUP {
    tag "$sample.id"
    label "process_medium_resources"
    label "tool_terminus"

    input:
        tuple val(sample), path(quant)

    output:
        path salmonOutput,      emit: salmon
        path output,            emit: terminus

    script:
        args = task.ext.args ?: ""
        output = "terminus/" + sample.id
        salmonOutput = "salmon." + sample.id

        """
        terminus group $args -d $quant -o terminus
        mv $quant $salmonOutput
        """

    stub:
        output = "terminus/" + sample.id
        salmonOutput = "salmon." + sample.id
        """
        mkdir -p $output
        touch $output/clusters.txt
        mv $quant $salmonOutput
        """
}
