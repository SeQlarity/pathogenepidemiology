include { DOWNLOAD_REFFASTA as DOWNLOAD_QUERYFASTA } from '../modules/local/download_reffasta/main'
include { DOWNLOAD_REFFASTA as DOWNLOAD_HOSTFASTA  } from '../modules/local/download_reffasta/main'
include { GUNZIP_REFFASTA                          } from '../modules/local/gunzip_reffasta/main'


// defaults are in nextflow.custom.config for now
//params.queryurl = null
//params.hosturl = null


workflow PREPARE_REFERENCES {

    take:
    queryurl
    hosturl

    main:
    ch_queryref = Channel.of(queryurl)
    ch_hostref  = Channel.of(hosturl)
    

    ch_query    = DOWNLOAD_QUERYFASTA(ch_queryref)
    hostfasta   = DOWNLOAD_HOSTFASTA(ch_hostref)
    queryfasta  = GUNZIP_REFFASTA(ch_query)

    emit:
    queryfasta
    hostfasta
}
