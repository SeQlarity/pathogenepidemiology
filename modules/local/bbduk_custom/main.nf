process BBDUK_CUSTOM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    // path changed to val to accept comma sep'd list like adapters,phix,path.fa
    path contaminants
    path adapters

    output:
    tuple val(meta), path('*.f*q*'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    tuple val(meta), path('*.discarded.f*q*'), emit: discarded
    tuple val(meta), path('*.stats.txt'), emit: stats
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_1.fastq.gz out2=${prefix}_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=${adapters},phix,${contaminants}" : "ref=${adapters},phix"
    def discarded = meta.single_end ? "outm=${prefix}.discarded.fastq.gz" : "outm1=${prefix}_1.discarded.fastq.gz outm2=${prefix}_2.discarded.fastq.gz"
    def stats_file = "stats=${prefix}.stats.txt"
    """
    bbduk.sh \\
        -Xmx${task.memory.toGiga()}g \\
        tossbrokenreads=t qin=33 \\
        $raw \\
        $trimmed \\
        $discarded \\
        $stats_file \\
        $contaminants_fa \\
        $args \\
        threads=$task.cpus \\
        &> ${prefix}.bbduk.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_command  = meta.single_end ? "echo '' | gzip > ${prefix}.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastq.gz ; echo '' | gzip > ${prefix}_2.fastq.gz"
    """
    touch ${prefix}.bbduk.log
    $output_command
    """
}
