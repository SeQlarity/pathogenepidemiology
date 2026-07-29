params.outdir = "${workflow.launchDir}/results/pipeline_results"
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
include { FASTQC             } from '../modules/nf-core/fastqc/main'
include { BBDUK_CUSTOM       } from '../modules/local/bbduk_custom/main'
include { MINIMAP2_INDEX     } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN     } from '../modules/nf-core/minimap2/align/main' 
include { MULTIQC            } from '../modules/nf-core/multiqc/main'
include { CLAIR3_CUSTOM      } from '../modules/local/clair3_custom/main' // CLAIR3 WILL BE PATCHED OUT FOR SOMETHING THAT HANDLES POLYPLOIDY BETTER
// use the CLAIR3 process when you want to use a locally-stored clair3 model 
   // Usage: tuple(meta, bam, bai, null, user_model, platform)
include { CLAIR3             } from '../modules/nf-core/clair3/main' // CLAIR3 WILL BE PATCHED OUT FOR SOMETHING THAT HANDLES POLYPLOIDY BETTER
include { BWAMEM3_INDEX      } from '../modules/nf-core/bwamem3/index/main'
include { BWAMEM3_MEM        } from '../modules/nf-core/bwamem3/mem/main'
include { SAMTOOLS_FAIDX     } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_MM2 } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_BM3 } from '../modules/nf-core/samtools/stats/main'


workflow {
  PREPARE_REFERENCES(
    params.queryurl,
    params.hosturl
  )

  ch_queryfasta = PREPARE_REFERENCES.out.queryfasta
  ch_hostfasta  = PREPARE_REFERENCES.out.hostfasta

  // Starting channels
  ch_samples = Channel
      .fromPath('testing/indir/samplesheet/samplesheet.csv')
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
      //.take(2) // for testing
  // Drop platform for FASTQC and BBDUK
  ch_reads = ch_samples.map { meta, reads, platform ->
      tuple(meta, reads)
  }

  ch_reffasta = ch_queryfasta
        .map { ref -> tuple([id: ref.baseName], ref) }
  ch_faidx_in = ch_reffasta
    .map { meta, fasta ->
        // Create a Path object for the expected .fai file
        def fai_path = file(fasta.toString() + ".fai")
        return tuple(meta, fasta, fai_path)
    }
  ch_reffai = SAMTOOLS_FAIDX(ch_faidx_in, false).fai

  // Produce QC reports per sample
  fastqc_out = FASTQC(ch_reads)
  

  ch_adapters = Channel.fromPath("${workflow.launchDir}/assets/adapters.fa")
  // Filter out host reads, adapters, phix
  hostrm_prepped = BBDUK_CUSTOM(ch_reads, ch_hostfasta.first(), ch_adapters.first())
  /*ch_mm2align_input = hostrm_prepped.reads.map { meta, reads -> 
        tuple(meta, reads, meta.single_end) 
    }*/
  
  ch_samples_l = hostrm_prepped.reads
    .join(ch_samples.map { meta, reads, platform -> tuple(meta, platform) }, by: 0)
    .filter { meta, reads, platform ->
        platform == 'OXFORD_NANOPORE'
    }
    .map { meta, reads, platform ->
        tuple(meta, reads)
    }
    //.take(5)

  ch_samples_s = hostrm_prepped.reads
    .join(ch_samples.map { meta, reads, platform -> tuple(meta, platform) }, by: 0)
    .filter { meta, reads, platform ->
        platform == 'ILLUMINA'
    }
    .map { meta, reads, platform ->
        tuple(meta, reads)
    }

  // Alignment of long reads
  minimap2_index = MINIMAP2_INDEX(ch_reffasta)
  aligned_l = MINIMAP2_ALIGN(
    ch_samples_l,                                                    
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

  // Alignment of short reads
  bwamem3_index = BWAMEM3_INDEX(ch_reffasta)
  aligned_s = BWAMEM3_MEM(
    ch_samples_s,                                                    
    bwamem3_index.index.first(),
    ch_reffasta.first(),
    true                                    // sort the bam for gatk                                                                    
    )

  // samtools stats input
  // takes tuple val(meta2), path(fasta), path(fai)
  statsrefs = ch_reffasta.join(ch_reffai, by: 0)
                          .map { meta, ref, fai_meta, fai -> tuple(meta, ref, fai) }
                          

  // Minimap2 stats
  ch_mm2stats_input = aligned_l.bam
    .join(aligned_l.index, by: 0)
    .map { meta, bam, bai -> tuple(meta, bam, bai) }


  mm2stats = SAMTOOLS_STATS_MM2(ch_mm2stats_input, statsrefs.first())

  // BWA-MEM3 stats
  ch_bm3stats_input = aligned_s.aligned
    .join(aligned_s.index, by: 0)
    .map { meta, bam, bai -> tuple(meta, bam, bai) }

  bm3stats = SAMTOOLS_STATS_BM3(ch_bm3stats_input, statsrefs.first())


  // multiqc 
  ch_multiqc_input = Channel.empty()
    .mix(
        fastqc_out.zip.collect { meta, files -> files },
        hostrm_prepped.stats.collect { meta, files -> files },
        hostrm_prepped.log.collect { meta, files -> files },
        hostrm_prepped.discarded.collect { meta, files -> files },
        mm2stats.stats.collect { meta, files -> files },
        bm3stats.stats.collect { meta, files -> files }
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
