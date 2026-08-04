
include { DOWNLOAD_REFFASTA as DOWNLOAD_QUERYFASTA } from '../../modules/local/download_reffasta/main'
include { DOWNLOAD_REFFASTA as DOWNLOAD_HOSTFASTA  } from '../../modules/local/download_reffasta/main'
include { GUNZIP_REFFASTA   as GUNZIP_QUERYFASTA   } from '../../modules/local/gunzip_reffasta/main'
include { SAMTOOLS_FAIDX                           } from '../../modules/nf-core/samtools/faidx/main'


workflow PREPARE_REFERENCES {

    take:
    queryurl
    hosturl

    main:
    ch_queryref = Channel.of(queryurl)
    ch_hostref  = Channel.of(hosturl)
    

    ch_query    = DOWNLOAD_QUERYFASTA(ch_queryref)
    hostfasta   = DOWNLOAD_HOSTFASTA(ch_hostref)
    queryfasta  = GUNZIP_QUERYFASTA(ch_query)

    queryfasta  = queryfasta.map { ref -> tuple([id: ref.baseName], ref) }
    
    ch_faidx_in = queryfasta
                    .map { meta, fasta ->
                        // Create a Path object for the expected .fai file
                        def fai_path = file(fasta.toString() + ".fai")
                        return tuple(meta, fasta, fai_path)
    }
    queryfai = SAMTOOLS_FAIDX(ch_faidx_in, false).fai

    emit:
    queryfasta
    queryfai
    hostfasta
}


  