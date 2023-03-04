import java.nio.file.Paths

def add_suffix(name, suffix) {
    name = name.toString()
    suffix != ""
                ? (
                    name.indexOf(".") != -1
                        ? name.substring(0, name.indexOf(".")) + suffix + name.substring(name.indexOf("."))
                        : name + task.ext.fileSuffix
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
        tuple val(sample), path(output), emit: trimmed
        path  report_json,               emit: json
        path  report_html,               emit: html
        path  "${task.ext.failedDir ?: 'removed-reads'}/*", optional: true, emit: orphaned
        path  "${task.ext.orphanedDir ?: 'orphaned-reads'}/*", optional: true, emit: filtered_out

    script:
        args   = task.ext.args ?: ""

        // Get output file names
        out1 = add_suffix(fastqs[0], task.ext.trimSuffix ?: "")
        out2 = add_suffix(fastqs[1] ?: "", task.ext.trimSuffix ?: "")
        output = [out1, out2]

        // Prepare input and output files
        input  = (sample.single == true ? "--in1 ${fastqs[0]}" : "--in1 ${fastqs[0]} --in2 ${fastqs[1]}")
        output = (sample.single == true ? "--out1 $out1" : "--out1 $out1 --out2 $out2")

        // Reports
        report_json = Paths.get((task.ext.reportDir ?: 'reports'), (task.ext.jsonDir ?: "json"), sample.id + ".json")
        report_html = Paths.get((task.ext.reportDir ?: 'reports'), (task.ext.htmlDir ?: "html"), sample.id + ".html")

        // Track discarded reads
        track_script = ""

        failed_file = Paths.get((task.ext.failedDir ?: 'removed-reads'), sample.id + (task.ext.removedPrefix ?: ".removed")  +".fastq.gz")
        orphan1 = Paths.get((task.ext.orphanedDir ?: 'orphaned-reads'), add_suffix(fastqs[0], task.ext.orphanSuffix ?: ""))
        orphan2 = Paths.get((task.ext.orphanedDir ?: 'orphaned-reads'), add_suffix(fastqs[1] ?: "", task.ext.orphanSuffix ?: ""))

        if (task.ext.trackDiscardedReads) {
            track_script = [
                "mkdir -p \$(dirname $failed_file) && touch $failed_file",
                "mkdir -p \$(dirname $orphan1) && touch $orphan1",
                fastqs[1] ? "mkdir -p \$(dirname $orphan2) && touch $orphan2" : ""
            ].join("\n").trim()
        }
        """
        fastp $input $output $track_script --thread $task.cpus --html $report_html --json $report_json $args
        """

    stub:
        // Get output file names
        out1 = add_suffix(fastqs[0], task.ext.trimSuffix ?: "")
        out2 = add_suffix(fastqs[1] ?: "", task.ext.trimSuffix ?: "")
        output = [out1, out2]

        report_json = Paths.get((task.ext.reportDir ?: 'reports'), (task.ext.jsonDir ?: "json"), sample.id + ".json")
        report_html = Paths.get((task.ext.reportDir ?: 'reports'), (task.ext.htmlDir ?: "html"), sample.id + ".html")

        // Track discarded reads
        track_script = ""

        failed_file = Paths.get((task.ext.failedDir ?: 'removed-reads'), sample.id + (task.ext.removedPrefix ?: ".removed")  +".fastq.gz")
        orphan1 = Paths.get((task.ext.orphanedDir ?: 'orphaned-reads'), add_suffix(fastqs[0], task.ext.orphanSuffix ?: ""))
        orphan2 = Paths.get((task.ext.orphanedDir ?: 'orphaned-reads'), add_suffix(fastqs[1] ?: "", task.ext.orphanSuffix ?: ""))

        if (task.ext.trackDiscardedReads) {
            track_script = [
                "mkdir -p \$(dirname $failed_file) && touch $failed_file",
                "mkdir -p \$(dirname $orphan1) && touch $orphan1",
                fastqs[1] ? "mkdir -p \$(dirname $orphan2) && touch $orphan2" : ""
            ].join("\n").trim()
        }
        """
        touch $out1
        ${fastqs[1] ? 'touch ' + out2 : ''}
        mkdir -p "\$(dirname $report_json)" && touch "$report_json"
        mkdir -p "\$(dirname $report_html)" && touch "$report_html"
        $track_script
        """
}
