nextflow.enable.dsl=2

process PHYLOGIC_SAMPLE {
    tag "$meta.sample_id"
    publishDir "${params.outdir}/samples", mode: 'copy', overwrite: false
    input:
    tuple val(meta), path(vcf), path(cn), path(purity)
    path drivers
    path peaks
    output:
    path '*.run'
    script:
    def encoded = groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    """
    python /opt/phylogic/workflow/run_nf_sample.py --meta '${encoded}' \
      --vcf '${vcf}' --cn '${cn}' --purity '${purity}' \
      --drivers '${drivers}' --peaks '${peaks}' \
      --iterations ${params.cluster_iterations} --seed ${params.seed}
    """
}

process COHORT_LEAGUE {
    publishDir "${params.outdir}", mode: 'copy', overwrite: false
    cpus 4
    memory '16 GB'
    time '24h'
    input:
    path runs
    output:
    path 'cohort'
    script:
    """
    python /opt/phylogic/workflow/run_nf_cohort.py \
      --expected ${params.expected_samples} --perms ${params.league_permutations} \
      --seasons ${params.league_seasons} --jobs ${task.cpus} --seed ${params.seed}
    """
}

workflow {
    if (!params.input || !params.outdir || !params.container || !params.driver_genes || !params.cn_peaks)
        error 'Require input, outdir, container, driver_genes and cn_peaks'
    rows = Channel.fromPath(params.input, checkIfExists: true).splitCsv(header: true).collect().map { items ->
        if (!items || items.size() != params.expected_samples as int)
            error "Expected ${params.expected_samples} samples; observed ${items.size()}"
        def ids = items.collect { it.sample_id }
        if (ids.unique(false).size() != items.size()) error 'Duplicate specimen IDs'
        items.each { r ->
            if (!(r.sample_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Unsafe specimen ID'
            ['vcf_path','cn_path','purity_path'].each { k -> if (!r[k]) error "Missing ${k}" }
            if (!(r.league_group ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Require explicit league_group; use EXCLUDED to omit from LeagueModel'
        }
        items
    }.flatMap { it }.map { r -> tuple(r, file(r.vcf_path, checkIfExists:true), file(r.cn_path, checkIfExists:true), file(r.purity_path, checkIfExists:true)) }
    PHYLOGIC_SAMPLE(rows, file(params.driver_genes, checkIfExists:true), file(params.cn_peaks, checkIfExists:true))
    COHORT_LEAGUE(PHYLOGIC_SAMPLE.out.collect())
}
