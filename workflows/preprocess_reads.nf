
include { FASTQC       } from '../modules/nf-core/fastqc/main'
include { BBDUK_CUSTOM } from '../modules/local/bbduk_custom/main'




workflow PREPROCESS_READS {

    take:
    samplesheet
    ch_hostfasta // from workflow PREPARE_REFERENCES
    adapters

    main:
    ch_samples = Channel.fromPath(samplesheet)
        .splitCsv(header: true, quote: '"')
        .map { row ->
            if (!row.fastq_1) throw new Exception("Missing fastq_1 for ${row.sample}")
            def meta = [
                id: row.run_accession,
                single_end: row.fastq_2 == '',
                instrument_platform: row.instrument_platform
                ]
            def reads = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
            tuple(meta, reads, row.instrument_platform)
            }
    // Drop platform for FASTQC and BBDUK
    ch_reads_raw = ch_samples.map { meta, reads, platform ->
        tuple(meta, reads)
    }

    

    // Produce QC reports per sample
    fastqc_out = FASTQC(ch_reads_raw)
    

    ch_adapters = Channel.fromPath(adapters)
    // Filter out host reads, adapters, phix
    ch_reads_clean = BBDUK_CUSTOM(ch_reads_raw, ch_hostfasta.first(), ch_adapters.first())

    emit:
    reads_clean   = ch_reads_clean.reads
    bbduk_stats   = ch_reads_clean.stats
    bbduk_logs    = ch_reads_clean.log
    bbduk_dropped = ch_reads_clean.discarded
    samples       = ch_samples
    fastqc        = fastqc_out.zip
}


