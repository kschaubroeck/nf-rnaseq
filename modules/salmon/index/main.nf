import java.nio.file.Files

process SALMON_INDEX {
    tag "$transcriptome_fasta"
    label "process_medium_resources"
    label "tool_salmon"

    input:
        path transcriptome_fasta
        path genome_fasta

    output:
        path 'salmon-index',      emit: index
        path "fasta/*",           emit: fasta
        path "decoy-names.txt",   optional: true, emit: decoys

    script:
        args = task.ext.args ?: ""

        // Determine if this is a gzipped file because if so then we need to extract the file
        gzipped      = Files.probeContentType(genome_fasta) == "application/gzip"
        genome_decoy = gzipped ? "<(gunzip -c $genome_fasta)" : $genome_fasta
        gentrome     = gzipped ? "gentrome.fa.gz" : "gentrome.fa"

        """
        mkdir fasta

        grep "^>" $genome_decoy | cut -d " " -f 1 > decoy-names.txt
        sed -i.bak -e 's/>//g' decoy-names.txt
        cat $transcriptome_fasta $genome_fasta > fasta/$gentrome

        mv -f $transcriptome_fasta fasta/$transcriptome_fasta
        mv -f $genome_fasta fasta/$genome_fasta

        salmon index --threads $task.cpus -t fasta/$gentrome -d decoy-names.txt --index salmon-index $args
        """

    stub:
        transcriptome_fasta = "fasta/transcriptome.fa.gz"
        genome_fasta        = "fasta/genome.fa.gz"
        gentrome            = "fasta/gentrome.fa.gz"
        """
        mkdir salmon-index
        mkdir fasta
        touch $transcriptome_fasta
        touch $genome_fasta
        touch $gentrome
        touch decoy-names.txt
        """
}
