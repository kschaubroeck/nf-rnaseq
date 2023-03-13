process TERMINUS_COLLAPSE {
    tag "$input"
    label "process_medium_resources"
    label "tool_terminus"

    input:
        path  quant,       stageAs: "salmon/*"
        path  terminus,    stageAs: "terminus/*" // Must match the default for terminus_dir

    output:
        path "$terminus_dir/*", includeInputs: true, emit: results

    script:
        args = task.ext.args ?: ""
        
        terminus_dir = task.ext.terminusDir ?: "terminus"
        terminus_code = (terminus_dir != "terminus") ? "mv terminus $terminus_dir" : ""

        input = (quant instanceof List ? quant : quant.toList()).join(' ')

        """
        terminus collapse $args -d $input -o terminus
        $terminus_code
        """

    stub:
        terminus_dir = task.ext.terminusDir ?: "terminus"
        terminus_code = (terminus_dir != "terminus") ? "mv terminus $terminus_dir" : ""

        input = (quant instanceof List ? quant : quant.toList()).join(' ')
        """
        mkdir -p terminus
        echo '$terminus $input' > terminus/collapsed.txt
        $terminus_code
        """
}
