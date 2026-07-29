process GUNZIP_REFFASTA {
    input:
    path ref_gz
    
    output:
    path "${ref_gz.baseName}"
    
    script:
    """
    gunzip -c ${ref_gz} > ${ref_gz.baseName}
    """
}