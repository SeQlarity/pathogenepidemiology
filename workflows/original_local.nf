nextflow.enable.dsl = 2


params.outdir = "${workflow.launchDir}/results"   
params.hostref = "${workflow.launchDir}/refs/Anopheles_gambiae.AgamP4.dna.toplevel.fa.gz"
params.query_genus = "Plasmodium"
params.query_species = "falciparum"
params.query_ref = "${workflow.launchDir}/refs/${params.query_genus}_${params.query_species}*.{fa,fasta,fa.gz,fasta.gz}"
params.threads = 4
params.use_gpu = false
// change --ont_qual to "hac" allowed
params.ont_qual = "sup"   
// allow user to fully or partially override clair3model value on cmd line, in case of older flowcell or "hac" basecalling
params.clair3model = params.clair3model ?: "r1041_e82_400bps_${params.ont_qual}_v520_with_mv"

include { FASTQC } from './modules/nf-core/fastqc/main'
include { BBDUK_CUSTOM } from './modules/local/bbduk_custom/main'
include { MINIMAP2_INDEX } from './modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from './modules/nf-core/minimap2/align/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { CLAIR3_CUSTOM } from './modules/local/clair3_custom/main'
// use the CLAIR3 process when you want to use a locally-stored clair3 model 
   // Usage: tuple(meta, bam, bai, null, user_model, platform)
include { CLAIR3 } from './modules/nf-core/clair3/main'



// wget link(s):
// reffasta: https://ftp.ebi.ac.uk/ensemblgenomes/pub/protists/release-62/fasta/plasmodium_falciparum/dna/Plasmodium_falciparum.GCA000002765v3.dna.toplevel.fa.gz

workflow {
  // Starting channels
  ch_samples = Channel
        .fromPath('indir/samplesheet/samplesheet.csv')
        .splitCsv(header: true, quote: '"')
        .map { row ->
            if (!row.fastq_1) throw new Exception("Missing fastq_1 for ${row.sample}")
            def meta = [id: row.run_accession, single_end: row.fastq_2 == '', instrument_platform: row.instrument_platform]
            def reads = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
            tuple(meta, reads, row.instrument_platform) 
        }
        //.take(2) // for testing
   
  ch_reads = ch_samples.map { meta, reads, platform -> 
      tuple(meta, reads) 
  }
  ch_hostref = Channel.of(params.hostref)
  ch_adapters = Channel.of("${workflow.launchDir}/adapters/adapters.fa")
  ch_reffasta = Channel.fromPath("${params.query_ref}")
    .map { ref -> tuple([id: ref.baseName], ref) }
  ch_reffai = ch_reffasta.map { meta, ref ->
    def fai = file(ref.toString() + ".fai")
    return tuple(meta, fai)
  } 

  // Produce QC reports per sample
  fastqc_out = FASTQC(ch_reads)
  
  // Filter out host reads, adapters, phix
  hostrm_prepped = BBDUK_CUSTOM(ch_reads, ch_hostref.first(), ch_adapters.first())
  ch_mm2align_input = hostrm_prepped.reads.map { meta, reads -> 
        tuple(meta, reads, meta.single_end) 
    }
  

  // Alignment of long reads
  minimap2_index = MINIMAP2_INDEX(ch_reffasta)
  aligned_l = MINIMAP2_ALIGN(
    hostrm_prepped.reads,                                                    
    minimap2_index.index.first(),      // reference as tuple
    true,                                                                    // bam_format
    'bai',                                                                   // bam_index_extension
    false,                                                                   // cigar_paf_format
    false                                                                    // cigar_bam
    )

  
  // Variant calling of long reads
  varcalls_l = CLAIR3_CUSTOM(
    aligned_l.bam.map { meta, bam ->
            def bai = file("${bam}.bai")
            def packaged_model = params.clair3model // must be null if using user_model
            def user_model = null // use process CLAIR3 if you are filling this input
            def platform = "ont"
            return tuple(meta, bam, bai, packaged_model, user_model, platform)
        },
    ch_reffasta.first(),
    ch_reffai.first()
  )


  // multiqc 
  ch_multiqc_input = Channel.empty()
    .mix(
        fastqc_out.zip.collect { meta, files -> files },
        hostrm_prepped.stats.collect { meta, files -> files },
        hostrm_prepped.log.collect { meta, files -> files },
        hostrm_prepped.discarded.collect { meta, files -> files }
        // Space for more channels
    )
    .flatten()
    .collect()
    .map { files ->
        def meta = [:]
        return tuple(meta, files, [], [], [], [])
    }
  MULTIQC(ch_multiqc_input)


}
