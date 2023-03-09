process TERMINUS_COLLAPSE {
    tag "$input"
    label "process_medium_resources"
    label "tool_terminus"

    input:
        path  quant,    stageAs: "salmon/*"
        path  terminus, stageAs: "terminus/*"

    output:
        path "terminus/*",  emit: results

    script:
        args = task.ext.args ?: ""
        input = (quant instanceof List ? quant : quant.toList()).join(' ')
        """
        terminus collapse $args -d $input -o terminus
        """

    stub:
        input = (quant instanceof List ? quant : quant.toList()).join(' ')
        """
        mkdir terminus
        echo '$terminus $input' > terminus/collapsed.txt
        """
}
