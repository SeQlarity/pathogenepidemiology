process DOWNLOAD_REFFASTA {
    storeDir "refs"
    
    input:
    val ref_url
    
    output:
    path "*.fa.gz"
    
    script:
    """
    wget -O ${ref_url.toString().split('/').last()} ${ref_url}
    """
}