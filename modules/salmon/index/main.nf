import java.nio.file.Files

process SALMON_INDEX {
    tag "$transcriptome_fasta"
    label "process_medium_resources"
    label "tool_salmon"

    input:
        path transcriptome_fasta
        path genome_fasta

    output:
        path "salmon-index",                emit: index
        path "*.{fasta,fasta.gz,fa,fa.gz}", includeInputs: true, emit: fasta
        path "decoy-names.txt",             optional: true, emit: decoys

    script:
        args = task.ext.args ?: ""
        index_name = task.ext.indexName ?: "salmon-index"

        // Determine if this is a gzipped file because if so then we need to extract the file
        gzipped      = Files.probeContentType(genome_fasta) == "application/gzip"
        genome_decoy = gzipped ? "<(gunzip -c $genome_fasta)" : $genome_fasta
        gentrome     = gzipped ? "gentrome.fa.gz" : "gentrome.fa"
        """
        grep "^>" $genome_decoy | cut -d " " -f 1 > decoy-names.txt
        sed -i.bak -e 's/>//g' decoy-names.txt
        cat $transcriptome_fasta $genome_fasta > $gentrome

        salmon index --threads $task.cpus -t $gentrome -d decoy-names.txt --index $index_name $args
        """
    stub:
        index_name = task.ext.indexName ?: "salmon-index"
        gentrome            = "gentrome.fa.gz"
        """
        mkdir $index_name
        touch $index_name/kmers.txt
        touch $gentrome
        touch decoy-names.txt
        """
}
