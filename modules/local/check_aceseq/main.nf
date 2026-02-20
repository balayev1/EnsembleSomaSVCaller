process CHECK_ACESEQ_DIR {
    executor 'local'
    
    input:
    val signal 

    exec:
    def aceseq_path = file(params.aceseq_outdir)
    
    println "\n------------------------------------------------"
    println "  Zero-shot CNV Calling Pipeline Finished. Verifying ACESeq Dir...    "
    println "------------------------------------------------"

    if (aceseq_path.exists()) {
        println "SUCCESS: Found pre-existing ACESeq directory at:"
        println "   ${aceseq_path}"
        println "------------------------------------------------\n"
    } else {
        println "ERROR: ACESeq directory NOT found at:"
        println "   ${aceseq_path}"
        println "------------------------------------------------\n"
        error "Missing ACESeq directory!"
    }
}