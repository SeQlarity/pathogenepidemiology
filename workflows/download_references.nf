// This will be patched out of first release, after which reference download will be fully automatic
// and done by ./prepare_references

include {DOWNLOAD_REFFASTA} from '../modules/local/download_reffasta/main'

params.queryurl = null
params.hosturl = null
//params.queryurl = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/protists/release-62/fasta/plasmodium_falciparum/dna/Plasmodium_falciparum.GCA000002765v3.dna.toplevel.fa.gz"
//params.hosturl = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/release-62/fasta/anopheles_gambiae/dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa.gz"

// TODO: expand to:
// workflow DOWNLOAD_REFERENCES
workflow {
    ch_refs = Channel.of(params.queryurl, params.hosturl)
    DOWNLOAD_REFFASTA(ch_refs)
}