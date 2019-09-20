// run as: nextflow smv_pipeline.nf -profile redwood --fasta /scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta --bams "/scratch/ucgd/lustre/UCGD_Datahub/Repository/AnalysisData/2018/A523/18-01-23_WashU-Yandell-CEPH_Sent/UCGD/GRCh38/Data/PolishedBams/*.cram" --project sfari

params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------

    smv_pipeline: a smoove workflow
    ============================

    smoove is available at: https://github.com/brentp/smoove

    Required arguments:
    -------------------

    --bams                Aligned sequences in .bam and/or .cram format.
                          Indexes (.bai/.crai) must be present.
    --fasta               Reference FASTA. Index (.fai) must exist in same
                          directory.
    Options:
    --------

    --outdir              Base results directory for output.
                          Default: '/.results'
    --project             File prefix for merged VCF files.
                          Default: 'sites'
    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// required arguments
params.fasta = false
if( !params.fasta ) { exit 1, "--fasta is not defined" }
params.bams = false
if( !params.bams ) { exit 1, "--bams is not defined" }

// variables
project = params.project ?: 'sites'
outdir = params.outdir ?: './results'
indexes = params.bams + ("${params.bams}".endsWith('.cram') ? '.crai' : '.bai')

log.info("\n")
log.info("Project            (--project)       : ${project}")
log.info("Reference fasta    (--fasta)         : ${params.fasta}")
log.info("Alignments         (--bams)          : ${params.bams}")
log.info("Indexes                              : ${indexes}")
log.info("Output             (--outdir)        : ${outdir}")
log.info("\n")

// instantiate files
fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")


// check file existence
if( !fasta.exists() ) { exit 1, "Missing reference fasta: ${fasta}" }
if( !faidx.exists() ) { exit 1, "Missing reference fasta index: ${faidx}" }


Channel
    .fromPath(params.bams, checkIfExists: true)
    .map { file -> tuple(file.baseName, file, file + ("${file}".endsWith('.cram') ? '.crai' : '.bai')) }
    .into { call_bams; genotype_bams }

Channel
    .fromPath(indexes, checkIfExists: true)
    .set { index_ch }

process smoove_call {
    publishDir path: "$outdir/smoove/called", mode: "copy", pattern: "*.vcf.gz*"
    publishDir path: "$outdir/logs", mode: "copy", pattern: "*-stats.txt"
    publishDir path: "$outdir/logs", mode: "copy", pattern: "*-smoove-call.log"

    input:
    set sample, file(bam), file(bai) from call_bams
    file fasta
    file faidx

    output:
    file("${sample}-smoove.genotyped.vcf.gz") into vcfs
    file("${sample}-smoove.genotyped.vcf.gz.csi") into idxs
    file("${sample}-stats.txt") into variant_counts
    file("${sample}-smoove-call.log") into sequence_counts

    script:
    """
    smoove call --genotype --name $sample --processes ${task.cpus} \
        --fasta $fasta \
        $bam 2> ${sample}-smoove-call.log
    bcftools stats ${sample}-smoove.genotyped.vcf.gz > ${sample}-stats.txt
    """
}

process smoove_merge {
    publishDir path: "$outdir/smoove/merged", mode: "copy"

    input:
    file vcf from vcfs.collect()
    file idx from idxs.collect()
    file fasta
    file faidx

    output:
    file("${project}.sites.vcf.gz") into sites

    script:
    """
    smoove merge --name $project --fasta $fasta $vcf
    """
}

process smoove_genotype {
    publishDir path: "$outdir/smoove/genotyped", mode: "copy"

    input:
    set sample, file(bam), file(bai) from genotype_bams
    file sites
    file fasta
    file faidx

    output:
    file("${sample}-smoove.genotyped.vcf.gz.csi") into genotyped_idxs
    file("${sample}-smoove.genotyped.vcf.gz") into genotyped_vcfs

    script:
    """
    wget -q https://raw.githubusercontent.com/samtools/samtools/develop/misc/seq_cache_populate.pl
    perl seq_cache_populate.pl -root \$(pwd)/cache $fasta 1> /dev/null 2> err || (cat err; exit 2)
    export REF_PATH=\$(pwd)/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=xx

    samtools quickcheck -v $bam
    smoove genotype --duphold --processes ${task.cpus} --removepr --outdir ./ --name ${sample} --fasta $fasta --vcf $sites $bam
    """
}

process smoove_square {
    publishDir path: "$outdir/smoove/squared", mode: "copy", pattern: "*.vcf.gz*"

    input:
    file vcf from genotyped_vcfs.collect()
    file idx from genotyped_idxs.collect()

    output:
    file("${project}.smoove.square.vcf.gz") into square_vcf
    file("${project}.smoove.square.vcf.gz.csi") into square_idx

    script:
    smoovepaste = "smoove paste --outdir ./ --name $project $vcf"
    if( vcf.collect().size() < 2 ) {
        paste = "cp $vcf ${project}.smoove.square.vcf.gz && cp $idx ${project}.smoove.square.vcf.gz.csi"
    }
    """
    $smoovepaste
    """
}
