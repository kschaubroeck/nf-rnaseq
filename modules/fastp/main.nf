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
        path  failed_file_path,                 optional: true, emit: failed
        tuple val(sample), path(orphan_paths),  optional: true, emit: orphans

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

        // Directory for failed and orphaned reads
        discard_dir = task.ext.discardDir ?: "discarded"

        // Suffixes for failed and orphaned reads
        orphanSuffix = task.ext.orphanSuffix ?: ".orphaned"
        failureSuffix = task.ext.removedSuffix ?: ".failed"

        // Names and paths for orphaned reads
        orphan1 = add_suffix(fastqs[0], orphanSuffix)
        orphan1_path = Paths.get(discard_dir, orphan1).toString()

        orphan2 = add_suffix(fastqs[1] ?: "", orphanSuffix)
        orphan2_path = Paths.get(discard_dir, orphan2).toString()

        orphan_paths = (sample.single == true ? [orphan1_path] : [orphan1_path, orphan2_path])

        // Failed (filtered out) reads paths
        failed_reads = add_suffix(sample.id, failureSuffix) + ".fastq.gz"
        failed_file_path = Paths.get(discard_dir, failed_reads).toString()

        // Get output file names
        out1 = add_suffix(fastqs[0], task.ext.trimSuffix ?: "")
        out2 = add_suffix(fastqs[1] ?: "", task.ext.trimSuffix ?: "")
        trimmed_names = [out1, out2]

        // Prepare input and output files CODE
        input  = "--in1 ${fastqs[0]}" + (sample.single ? "" : " --in2 ${fastqs[1]}")
        output = "--out1 $out1" + (sample.single? "" : " --out2 $out2")

        // Assemble the code for trakcing discarded reads 
        track_discard = task.ext.trackDiscardedReads ?: false  
        discard = "--failed_out $failed_reads --unpaired1 $orphan1" + (sample.single ? "" : " --unpaired2 $orphan2")

        // Check if we need to track the discards. If not, remove the code to do it
        if (track_discard == false) {
            discard = ""
        }
        """
        fastp $input $output $discard --thread $task.cpus --html $html_report_file --json $json_report_file $args

        mkdir $report_dir
        mkdir $json_dir
        mkdir $html_dir
        mkdir $discard_dir

        # Move items to their correct locations ONLY if they exists
        [ -f $json_report_file ] && mv $json_report_file $json_report
        [ -f $html_report_file ] && mv $html_report_file $html_report
        [ -f $orphan1 ] && mv $orphan1 $orphan1_path
        [ -f $orphan2 ] && mv $orphan2 $orphan2_path
        [ -f $failed_reads ] && mv $failed_reads $failed_file_path
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

        // Directory for failed and orphaned reads
        discard_dir = task.ext.discardDir ?: "discarded"

        // Suffixes for failed and orphaned reads
        orphanSuffix = task.ext.orphanSuffix ?: ".orphaned"
        failureSuffix = task.ext.removedSuffix ?: ".failed"

        // Names and paths for orphaned reads
        orphan1 = add_suffix(fastqs[0], orphanSuffix)
        orphan1_path = Paths.get(discard_dir, orphan1).toString()

        orphan2 = add_suffix(fastqs[1] ?: "", orphanSuffix)
        orphan2_path = Paths.get(discard_dir, orphan2).toString()

        orphan_paths = (sample.single == true ? [orphan1_path] : [orphan1_path, orphan2_path])

        // Failed (filtered out) reads paths
        failed_reads = add_suffix(sample.id, failureSuffix) + ".fastq.gz"
        failed_file_path = Paths.get(discard_dir, failed_reads).toString()

        // Get output file names
        out1 = add_suffix(fastqs[0], task.ext.trimSuffix ?: "")
        out2 = add_suffix(fastqs[1] ?: "", task.ext.trimSuffix ?: "")
        trimmed_names = [out1, out2]

        // Assemble the code for trakcing discarded reads 
        track_discard = task.ext.trackDiscardedReads ?: false  

        """
        mkdir report_dir
        mkdir json_dir
        mkdir html_dir
        mkdir discard_dir

        touch $out1
        ${fastqs[1] ? 'touch ' + out2 : ''}

        touch $html_report
        touch $json_report

        ${track_discard ? 'touch ' + failed_file_path : ''}
        ${track_discard ? 'touch ' + orphan1 : ''}
        ${(track_discard) ? (fastqs[1] ? 'touch ' + orphan2 : '') : ''}
        """
}
