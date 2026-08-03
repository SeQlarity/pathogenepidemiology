
//params.query_genus = "Plasmodium"     //not used for now
//params.query_species = "falciparum"   //not used for now
//params.query_ref = "${params.outdir}/download/${params.query_genus}_${params.query_species}*.{fa,fasta,fa.gz,fasta.gz}" //not used for now
params.threads = 4


// CLAIR3 WILL BE PATCHED OUT FOR SOMETHING THAT HANDLES POLYPLOIDY BETTER

// change --ont_qual to "hac" allowed
params.ont_qual = "sup"   
// allow user to fully or partially override clair3model value on cmd line, in case of older flowcell or "hac" basecalling
params.clair3model = params.clair3model ?: "r1041_e82_400bps_${params.ont_qual}_v520_with_mv"

// CLAIR3 WILL BE PATCHED OUT FOR SOMETHING THAT HANDLES POLYPLOIDY BETTER


include { PREPARE_REFERENCES } from './prepare_references.nf'
include { PREPROCESS_READS   } from './preprocess_reads.nf'
include { MINIMAP2_INDEX     } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN     } from '../modules/nf-core/minimap2/align/main' 
include { MULTIQC            } from '../modules/nf-core/multiqc/main'
include { CLAIR3_CUSTOM      } from '../modules/local/clair3_custom/main' // CLAIR3 WILL BE PATCHED OUT FOR SOMETHING THAT HANDLES POLYPLOIDY BETTER
// use the CLAIR3 process when you want to use a locally-stored clair3 model 
   // Usage: tuple(meta, bam, bai, null, user_model, platform)
include { CLAIR3             } from '../modules/nf-core/clair3/main' // CLAIR3 WILL BE PATCHED OUT FOR SOMETHING THAT HANDLES POLYPLOIDY BETTER
include { BWAMEM3_INDEX      } from '../modules/nf-core/bwamem3/index/main'
include { BWAMEM3_MEM        } from '../modules/nf-core/bwamem3/mem/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_MM2 } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_BM3 } from '../modules/nf-core/samtools/stats/main'


workflow {
  // DOWNLOAD AND INDEX REFERENCE FASTA FILES
  PREPARE_REFERENCES(
    params.queryurl,
    params.hosturl
  )

  ch_queryfasta = PREPARE_REFERENCES.out.queryfasta
  ch_hostfasta  = PREPARE_REFERENCES.out.hostfasta
  ch_reffai     = PREPARE_REFERENCES.out.queryfai


  // RUN FASTQC ON READS AND USE BBDUK TO REMOVE ADAPTER AND BIO. CONTAMINANT SEQS
  PREPROCESS_READS(
    params.samplesheet,
    ch_hostfasta,
    params.adapters
  )

  ch_reads_clean    = PREPROCESS_READS.out.reads_clean
  ch_bbduk_stats    = PREPROCESS_READS.out.bbduk_stats
  ch_bbduk_logs     = PREPROCESS_READS.out.bbduk_logs
  ch_bbduk_dropped  = PREPROCESS_READS.out.bbduk_dropped
  ch_samples        = PREPROCESS_READS.out.samples
  ch_fastqc_out     = PREPROCESS_READS.out.fastqc



  // LONG READ TRACK - ALIGNMENT + VAR CALLING

  // Filter short reads out
  ch_samples_l = ch_reads_clean
    .join(ch_samples.map { meta, reads, platform -> tuple(meta, platform) }, by: 0)
    .filter { meta, reads, platform ->
        platform == 'OXFORD_NANOPORE'
    }
    .map { meta, reads, platform ->
        tuple(meta, reads)
    }
    

  // Alignment of long reads
  minimap2_index = MINIMAP2_INDEX(ch_queryfasta)
  //minimap2_index = MINIMAP2_INDEX(ch_queryfasta)
  aligned_l = MINIMAP2_ALIGN(
    ch_samples_l,                                                    
    minimap2_index.index.first(),      // reference as tuple
    true,                                                                    // bam_format
    'bai',                                                                   // bam_index_extension
    false,                                                                   // cigar_paf_format
    false                                                                    // cigar_bam
    )

  // Minimap2 stats
  ch_mm2stats_input = aligned_l.bam
    .join(aligned_l.index, by: 0)
    .map { meta, bam, bai -> tuple(meta, bam, bai) }

  // Variant calling of long reads
  varcalls_l = CLAIR3_CUSTOM(
    aligned_l.bam.map { meta, bam ->
            def bai = file("${bam}.bai")
            def packaged_model = params.clair3model // must be null if using user_model
            def user_model = null // use process CLAIR3 if you are filling this input
            def platform = "ont"
            return tuple(meta, bam, bai, packaged_model, user_model, platform)
        },
    ch_queryfasta.first(),
    //ch_queryfasta.first(),
    ch_reffai.first()
  )





  // SHORT READ TRACK - ALIGNMENT + VAR CALLING

  // Filter long reads out
  ch_samples_s = ch_reads_clean
    .join(ch_samples.map { meta, reads, platform -> tuple(meta, platform) }, by: 0)
    .filter { meta, reads, platform ->
        platform == 'ILLUMINA'
    }
    .map { meta, reads, platform ->
        tuple(meta, reads)
    }

  // Alignment of short reads
  bwamem3_index = BWAMEM3_INDEX(ch_queryfasta)
  aligned_s = BWAMEM3_MEM(
    ch_samples_s,                                                    
    bwamem3_index.index.first(),
    ch_queryfasta.first(),
    //ch_queryfasta.first(),
    true                                    // sort the bam for gatk                                                                    
    )

  // BWA-MEM3 stats
  ch_bm3stats_input = aligned_s.aligned
    .join(aligned_s.index, by: 0)
    .map { meta, bam, bai -> tuple(meta, bam, bai) }


  // TODO: CREATE varcalls_s CHANNEL WITH GATK

  // TODO: MERGE varcalls_l AND varcalls_s CHANNELS




  // POST PROCESSING (eg final ops for multiqc report)

  // samtools stats input
  // takes tuple val(meta2), path(fasta), path(fai)
  
  statsrefs = ch_queryfasta.join(ch_reffai, by: 0)
                          .map { meta, ref, fai_meta, fai -> tuple(meta, ref, fai) }
                          
  ch_mm2stats = SAMTOOLS_STATS_MM2(ch_mm2stats_input, statsrefs.first())
  ch_bm3stats = SAMTOOLS_STATS_BM3(ch_bm3stats_input, statsrefs.first())

  // multiqc 
  ch_multiqc_input = Channel.empty()
    .mix(
        ch_fastqc_out.collect { meta, files -> files },
        ch_bbduk_stats.collect { meta, files -> files },
        ch_bbduk_logs.collect { meta, files -> files },
        ch_bbduk_dropped.collect { meta, files -> files },
        ch_mm2stats.stats.collect { meta, files -> files },
        ch_bm3stats.stats.collect { meta, files -> files }
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
