import java.nio.file.Paths

def add_suffix(name, suffix) {
    name = name.toString()
    suffix != ""
                ? (
                    name.indexOf(".") != -1
                        ? name.substring(0, name.indexOf(".")) + suffix + name.substring(name.indexOf("."))
                        : name + suffix
                )
                : name
}

process FASTP {
    tag "$sample.id"
    label "process_medium_resources"
    label "tool_fastp"

    input:
        tuple val(sample), path(fastqs)

    output:
        tuple val(sample), path(trimmed_names), emit: trimmed
        path  json_report,                      emit: json
        path  html_report,                      emit: html

    script:
        args   = task.ext.args ?: ""

        // Reports
        report_dir = task.ext.reportDir ?: "reports"

        json_dir = Paths.get(report_dir, task.ext.jsonDir ?: "json").toString()
        json_report_file = sample.id + ".json"
        json_report = Paths.get(json_dir, json_report_file).toString()

        html_dir = Paths.get(report_dir, task.ext.htmlDir ?: "html").toString()
        html_report_file = sample.id + ".html"
        html_report = Paths.get(html_dir, html_report_file).toString()

        // Get output file names
        out1 = add_suffix(fastqs[0], task.ext.trimSuffix ?: "")
        out2 = add_suffix(fastqs[1] ?: "", task.ext.trimSuffix ?: "")
        trimmed_names = [out1, out2]

        // Prepare input and output files CODE
        input  = "--in1 ${fastqs[0]}" + (sample.single ? "" : " --in2 ${fastqs[1]}")
        output = "--out1 $out1" + (sample.single? "" : " --out2 $out2")

        """
        fastp $input $output --thread $task.cpus --html $html_report_file --json $json_report_file $args

        mkdir $report_dir
        mkdir $json_dir
        mkdir $html_dir

        mv $json_report_file $json_report
        mv $html_report_file $html_report
        """

    stub:
        // Reports
        report_dir = task.ext.reportDir ?: "reports"

        json_dir = Paths.get(report_dir, task.ext.jsonDir ?: "json").toString()
        json_report_file = sample.id + ".json"
        json_report = Paths.get(json_dir, json_report_file).toString()

        html_dir = Paths.get(report_dir, task.ext.htmlDir ?: "html").toString()
        html_report_file = sample.id + ".html"
        html_report = Paths.get(html_dir, html_report_file).toString()

        // Get output file names
        out1 = add_suffix(fastqs[0], task.ext.trimSuffix ?: "")
        out2 = add_suffix(fastqs[1] ?: "", task.ext.trimSuffix ?: "")
        trimmed_names = [out1, out2]
        """
        mkdir $report_dir
        mkdir $json_dir
        mkdir $html_dir

        touch $out1
        ${fastqs[1] ? 'touch ' + out2 : ''}

        touch $html_report
        touch $json_report
        """
}
