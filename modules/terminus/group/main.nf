process TERMINUS_GROUP {
    tag "$sample.id"
    label "process_medium_resources"
    label "tool_terminus"

    input:
        tuple val(sample), path(quant)

    output:
        path salmon_output,     emit: salmon
        path "terminus/*",      emit: terminus

    script:
        args = task.ext.args ?: ""
        salmon_output = quant
        """
        terminus group $args -d $salmon_output -o terminus
        """

    stub:
        args = task.ext.args ?: ""
        salmon_output = quant
        """
        mkdir -p terminus
        touch terminus/${sample.id}
        """
}
