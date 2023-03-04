process TERMINUS_COLLAPSE {
    tag "$input"
    label "process_medium_resources"
    label "tool_terminus"

    input:
        path  quant
        path  terminus

    output:
        path "collapsed/*", emit: results
        path terminus,      emit: input

    script:
        args = task.ext.args ?: ""
        input = (quant instanceof List ? quant : quant.toList()).join(' ')
        """
        terminus collapse $args -d $input -o collapsed
        """

    stub:
        input = (quant instanceof List ? quant : quant.toList()).join(' ')
        """
        mkdir -p collapsed
        echo '$terminus $input' > collapsed/collapsed.txt
        touch collapsed/seconditem.txt
        """
}
