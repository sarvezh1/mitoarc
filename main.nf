nextflow.enable.dsl=2

include { INPUT_CHECK } from './modules/input_check'
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { PREPARE_REFERENCE } from './modules/prepare_reference'
include { PREPARE_SAMTOOLS_INDEX } from './modules/prepare_samtools_index'
include { PREPARE_DOWNSTREAM_REFERENCE } from './modules/prepare_downstream_reference'
include { BWA_MEM2 } from './modules/bwa_mem2'
include { SORT_INDEX } from './modules/sort_index'
include { FLAGSTAT as WHOLE_FLAGSTAT } from './modules/flagstat'
include { FLAGSTAT as INPUT_FLAGSTAT } from './modules/flagstat'
include { EXTRACT_MTDNA } from './modules/extract_mtdna'
include { FLAGSTAT as MTDNA_FLAGSTAT } from './modules/flagstat'
include { SELECT_AUTOSOMAL_CONTIGS } from './modules/select_autosomal_contigs'
include { MTDNA_COVERAGE } from './modules/mtdna_coverage'
include { AUTOSOMAL_COVERAGE } from './modules/autosomal_coverage'
include { COMBINE_MTDNA_METRICS } from './modules/combine_mtdna_metrics'
include { PREPARE_MT_REFERENCES } from './modules/prepare_mt_references'
include { PREPARE_MT_READS } from './modules/prepare_mt_reads'
include { REALIGN_MT as REALIGN_MT_STANDARD; REALIGN_MT as REALIGN_MT_SHIFTED } from './modules/realign_mt'
include { SORT_REALIGN_MT as SORT_REALIGN_MT_STANDARD; SORT_REALIGN_MT as SORT_REALIGN_MT_SHIFTED } from './modules/sort_realign_mt'
include { REALIGN_MT_FLAGSTAT as REALIGN_MT_STANDARD_FLAGSTAT; REALIGN_MT_FLAGSTAT as REALIGN_MT_SHIFTED_FLAGSTAT } from './modules/realign_mt_flagstat'
include { PREPARE_GATK_MT_REFERENCES } from './modules/prepare_gatk_mt_references'
include { GATK_MUTECT2 as GATK_MUTECT2_STANDARD; GATK_MUTECT2 as GATK_MUTECT2_SHIFTED } from './modules/gatk_mutect2'
include { LIFTOVER_MT_VCF } from './modules/liftover_mt_vcf'
include { COMBINE_GATK_VCFS } from './modules/combine_gatk_vcfs'
include { FILTER_GATK_FINAL } from './modules/filter_gatk_final'
include { VALIDATE_GATK_FIXTURE } from './modules/validate_gatk_fixture'
include { MUTSERVE2_CALL } from './modules/mutserve2_call'
include { NORMALIZE_MUTSERVE2 } from './modules/normalize_mutserve2'
include { VALIDATE_MUTSERVE2_FIXTURE } from './modules/validate_mutserve2_fixture'
include { COMPARE_CALLERS } from './modules/compare_callers'
include { NUMT_MAPPING_FILTER } from './modules/numt_mapping_filter'
include { PREPARE_MT_READS as PREPARE_MT_READS_B } from './modules/prepare_mt_reads'
include { REALIGN_MT as REALIGN_MT_STANDARD_B; REALIGN_MT as REALIGN_MT_SHIFTED_B } from './modules/realign_mt'
include { SORT_REALIGN_MT as SORT_REALIGN_MT_STANDARD_B; SORT_REALIGN_MT as SORT_REALIGN_MT_SHIFTED_B } from './modules/sort_realign_mt'
include { GATK_MUTECT2 as GATK_MUTECT2_STANDARD_B; GATK_MUTECT2 as GATK_MUTECT2_SHIFTED_B } from './modules/gatk_mutect2'
include { LIFTOVER_MT_VCF as LIFTOVER_MT_VCF_B } from './modules/liftover_mt_vcf'
include { COMBINE_GATK_VCFS as COMBINE_GATK_VCFS_B } from './modules/combine_gatk_vcfs'
include { FILTER_GATK_FINAL as FILTER_GATK_FINAL_B } from './modules/filter_gatk_final'
include { MUTSERVE2_CALL as MUTSERVE2_CALL_B } from './modules/mutserve2_call'
include { NORMALIZE_MUTSERVE2 as NORMALIZE_MUTSERVE2_B } from './modules/normalize_mutserve2'
include { NUMT_FIXTURE_PREFLIGHT } from './modules/numt_fixture_preflight'
include { NUMT_EVIDENCE } from './modules/numt_evidence'
include { GATK_NUMT_AUDIT } from './modules/gatk_numt_audit'
include { HAPLOGREP_CLASSIFY as HAPLOGREP_GATK; HAPLOGREP_CLASSIFY as HAPLOGREP_MUTSERVE2; HAPLOGREP_CLASSIFY as HAPLOGREP_FUNCTIONAL } from './modules/haplogrep_classify'
include { COMBINE_HAPLOGROUPS } from './modules/combine_haplogroups'
include { HAPLOCHECK_CONTAMINATION as HAPLOCHECK_GATK; HAPLOCHECK_CONTAMINATION as HAPLOCHECK_MUTSERVE2 } from './modules/haplocheck_contamination'
include { COMBINE_CONTAMINATION } from './modules/combine_contamination'
include { VALIDATE_HAPLOGREP_FUNCTIONAL } from './modules/validate_haplogrep_functional'
include { ANNOTATE_MTDNA as ANNOTATE_GATK; ANNOTATE_MTDNA as ANNOTATE_MUTSERVE2 } from './modules/annotate_mtdna'
include { CONSENSUS_MTDNA as CONSENSUS_GATK; CONSENSUS_MTDNA as CONSENSUS_MUTSERVE2 } from './modules/consensus_mtdna'
include { COMBINE_ANNOTATION_SUMMARIES } from './modules/combine_annotation_summaries'
include { VALIDATE_ANNOTATION_FUNCTIONAL } from './modules/validate_annotation_functional'
include { REPORTING_DEPTH } from './modules/reporting_depth'
include { PREPARE_REPORTING_DATA } from './modules/prepare_reporting_data'
include { GENERATE_REPORT_PLOTS } from './modules/generate_report_plots'
include { RENDER_MTDNA_REPORT } from './modules/render_mtdna_report'
include { VALIDATE_REPORTING } from './modules/validate_reporting'
include { MULTIQC_REPORT } from './modules/multiqc_report'
include { VALIDATE_MULTIQC_REPORT } from './modules/validate_multiqc_report'
include { REPORTING_FUNCTIONAL_FIXTURE } from './modules/reporting_functional_fixture'
include { VALIDATE_PREALIGNED } from './modules/validate_pre_aligned'
include { ALIGNMENT_METADATA } from './modules/alignment_metadata'
include { FINALIZE_ALIGNMENT_CONTRACT } from './modules/finalize_alignment_contract'

params.input = null
params.reference = null
params.mt_contig = 'chrM'
params.autosomal_contigs = null
params.coverage_mapq = 20
params.shift_offset = 8000
params.fixture_truth = null
params.mutserve2_heteroplasmy_threshold = 0.01
params.numt_mapq = 20
params.numt_strategy = 'none'
params.numt_fixture_truth = null
params.numt_fixture_design = null
params.haplogrep_tree = 'phylotree-rcrs@17.3'
params.haplogrep_het_level = 0.9
params.haplogroup_functional_fixture = null
params.haplogroup_functional_expected = null
params.homoplasmy_threshold = 0.95
params.consensus_heteroplasmy_mode = 'reference'
params.annotation_functional_fixture = null
params.annotation_functional_expected = null

workflow {

    log.info """
    ============================================
                     MitoArc
    Human mitochondrial genome analysis
    ============================================
    Nextflow version : ${workflow.nextflow.version}
    Run name         : ${workflow.runName}
    Project directory: ${projectDir}
    ============================================
    """

    if (!params.input) {
        error "Missing required parameter: --input"
    }
    if (!params.reference) {
        error "Missing required parameter: --reference"
    }
    if (!params.mt_contig) {
        error "Missing required parameter: --mt_contig"
    }
    if (!(params.coverage_mapq.toString() ==~ /[0-9]+/)) {
        error "Invalid --coverage_mapq value '${params.coverage_mapq}'. Use a non-negative integer."
    }
    if (!(params.shift_offset.toString() ==~ /[0-9]+/) || params.shift_offset.toString().toBigInteger() <= 0) {
        error "Invalid --shift_offset value '${params.shift_offset}'. Use an integer greater than zero."
    }
    if (params.autosomal_contigs && !(params.autosomal_contigs.toString() ==~ /[A-Za-z0-9_.-]+(,[A-Za-z0-9_.-]+)*/)) {
        error "Invalid --autosomal_contigs value. Use a comma-separated list of contig names."
    }
    if (params.fixture_truth && !file(params.fixture_truth).exists()) {
        error "Fixture truth file does not exist: ${params.fixture_truth}"
    }
    if (!(params.mutserve2_heteroplasmy_threshold.toString() ==~ /(?:0(?:\.\d+)?|1(?:\.0+)?)$/) || params.mutserve2_heteroplasmy_threshold.toBigDecimal() <= 0 || params.mutserve2_heteroplasmy_threshold.toBigDecimal() > 1) {
        error "Invalid --mutserve2_heteroplasmy_threshold value '${params.mutserve2_heteroplasmy_threshold}'. Use a number greater than 0 and at most 1."
    }
    if (!(params.numt_mapq.toString() ==~ /[0-9]+/) || params.numt_mapq.toString().toBigInteger() > 255) {
        error "Invalid --numt_mapq value '${params.numt_mapq}'. Use an integer from 0 through 255."
    }
    if (!(params.numt_strategy in ['none', 'comparison'])) {
        error "Invalid --numt_strategy value '${params.numt_strategy}'. Use none or comparison."
    }
    if (params.numt_fixture_truth && !file(params.numt_fixture_truth).exists()) {
        error "NUMT fixture truth file does not exist: ${params.numt_fixture_truth}"
    }
    if (params.numt_fixture_design && !file(params.numt_fixture_design).exists()) {
        error "NUMT fixture design file does not exist: ${params.numt_fixture_design}"
    }
    if (params.haplogrep_tree != 'phylotree-rcrs@17.3') {
        error "This pinned container provides --haplogrep_tree phylotree-rcrs@17.3 only."
    }
    if (!(params.haplogrep_het_level.toString() ==~ /(?:0(?:\.\d+)?|1(?:\.0+)?)$/)) {
        error "Invalid --haplogrep_het_level value '${params.haplogrep_het_level}'. Use a number from 0 through 1."
    }
    if ((params.haplogroup_functional_fixture && !params.haplogroup_functional_expected) || (!params.haplogroup_functional_fixture && params.haplogroup_functional_expected)) {
        error "Provide both --haplogroup_functional_fixture and --haplogroup_functional_expected, or neither."
    }
    if (params.haplogroup_functional_fixture && !file(params.haplogroup_functional_fixture).exists()) {
        error "Haplogroup functional fixture does not exist: ${params.haplogroup_functional_fixture}"
    }
    if (params.haplogroup_functional_expected && !file(params.haplogroup_functional_expected).exists()) {
        error "Haplogroup functional expectation does not exist: ${params.haplogroup_functional_expected}"
    }
    if (!(params.homoplasmy_threshold.toString() ==~ /(?:0(?:\.\d+)?|1(?:\.0+)?)$/) || params.homoplasmy_threshold.toBigDecimal() <= 0) {
        error "Invalid --homoplasmy_threshold value '${params.homoplasmy_threshold}'. Use a number greater than 0 and at most 1."
    }
    if (!(params.consensus_heteroplasmy_mode in ['reference', 'iupac'])) {
        error "Invalid --consensus_heteroplasmy_mode value '${params.consensus_heteroplasmy_mode}'. Use reference or iupac."
    }
    if ((params.annotation_functional_fixture && !params.annotation_functional_expected) || (!params.annotation_functional_fixture && params.annotation_functional_expected)) {
        error "Provide both --annotation_functional_fixture and --annotation_functional_expected, or neither."
    }
    if (params.annotation_functional_fixture && !file(params.annotation_functional_fixture).exists()) {
        error "Annotation functional fixture does not exist: ${params.annotation_functional_fixture}"
    }
    if (params.annotation_functional_expected && !file(params.annotation_functional_expected).exists()) {
        error "Annotation functional expectation does not exist: ${params.annotation_functional_expected}"
    }

    samplesheet = file(params.input)
    if (!samplesheet.exists()) {
        error "Samplesheet does not exist: ${params.input}"
    }

    reference = file(params.reference)
    if (!reference.exists()) {
        error "Reference FASTA does not exist: ${params.reference}"
    }

    header = samplesheet.readLines().find { it.trim() }
    if (!header) {
        error "Samplesheet is empty: ${params.input}"
    }

    columns = header.split(',', -1)*.trim()
    schemas = [
        ['sample', 'fastq_1', 'fastq_2'],
        ['sample', 'bam'],
        ['sample', 'cram']
    ]
    if (!schemas.any { it == columns }) {
        error "Invalid samplesheet schema. Use exactly sample,fastq_1,fastq_2 or sample,bam or sample,cram; mixed schemas are not supported. Found: ${columns.join(',')}"
    }
    input_mode = columns == ['sample', 'bam'] ? 'BAM' : (columns == ['sample', 'cram'] ? 'CRAM' : 'FASTQ')

    def seen_samples = [] as Set
    samples_ch = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def sample = row.sample?.toString()?.trim()
            def fq1_value = row.fastq_1?.toString()?.trim()
            def fq2_value = row.fastq_2?.toString()?.trim()
            def alignment_value = input_mode == 'BAM' ? row.bam?.toString()?.trim() : row.cram?.toString()?.trim()

            if (!sample) {
                error "Samplesheet contains an empty sample ID."
            }
            if (seen_samples.contains(sample)) {
                error "Samplesheet contains duplicate sample ID '${sample}'."
            }
            seen_samples.add(sample)

            if (input_mode != 'FASTQ') {
                if (!alignment_value) error "Sample '${sample}': ${input_mode.toLowerCase()} cannot be empty."
                def expected_suffix = input_mode == 'BAM' ? /(?i).*\.bam$/ : /(?i).*\.cram$/
                if (!(alignment_value ==~ expected_suffix)) error "Sample '${sample}': ${input_mode} samplesheet requires a .${input_mode.toLowerCase()} file (found '${alignment_value}')."
                def source_path = alignment_value.startsWith('/') ? file(alignment_value) : file(samplesheet.parent.resolve(alignment_value).normalize())
                if (!source_path.exists()) source_path = file(projectDir.resolve(alignment_value).normalize())
                if (!source_path.exists()) error "Sample '${sample}': ${input_mode} file does not exist: ${alignment_value}"
                def index_suffixes = input_mode == 'BAM' ? ['.bai', '.csi'] : ['.crai']
                def index_candidates = index_suffixes.collect { suffix -> file("${source_path}${suffix}") }
                if (!alignment_value.startsWith('/')) index_candidates += index_suffixes.collect { suffix -> file(samplesheet.parent.resolve("${alignment_value}${suffix}").normalize().toString()) }
                def input_index = index_candidates.find { it.exists() } ?: file("${projectDir}/assets/no_index")
                return tuple(sample, input_mode, source_path, reference, input_index)
            }
            if (!fq1_value || !fq2_value) error "Sample '${sample}': fastq_1 and fastq_2 cannot be empty."
            if (!(fq1_value ==~ /.*\.(fastq|fq)(\.gz)?/)) {
                error "Sample '${sample}': fastq_1 must end in .fastq, .fq, .fastq.gz, or .fq.gz (found '${fq1_value}')."
            }
            if (!(fq2_value ==~ /.*\.(fastq|fq)(\.gz)?/)) {
                error "Sample '${sample}': fastq_2 must end in .fastq, .fq, .fastq.gz, or .fq.gz (found '${fq2_value}')."
            }

            def fq1_path = fq1_value.startsWith('/') ? file(fq1_value) : file(projectDir.resolve(fq1_value).normalize())
            def fq2_path = fq2_value.startsWith('/') ? file(fq2_value) : file(projectDir.resolve(fq2_value).normalize())
            def fq1 = fq1_path.exists() ? fq1_path : file(samplesheet.parent.resolve(fq1_value).normalize())
            def fq2 = fq2_path.exists() ? fq2_path : file(samplesheet.parent.resolve(fq2_value).normalize())
            if (!fq1.exists()) {
                error "Sample '${sample}': FASTQ file does not exist: ${fq1_value}"
            }
            if (!fq2.exists()) {
                error "Sample '${sample}': FASTQ file does not exist: ${fq2_value}"
            }
            if (fq1.toAbsolutePath().normalize() == fq2.toAbsolutePath().normalize()) {
                error "Sample '${sample}': fastq_1 and fastq_2 must point to different files."
            }

            tuple(sample, fq1, fq2)
        }

    if (input_mode == 'FASTQ') {
        raw_reads = INPUT_CHECK(samples_ch)
        fastqc_reports = FASTQC(raw_reads)
        fastp_result = FASTP(fastqc_reports)
        alignment_reference_bundle = PREPARE_REFERENCE(Channel.value(reference))
        downstream_reference_bundle = PREPARE_SAMTOOLS_INDEX(alignment_reference_bundle)
        bwa_inputs = fastp_result.cleaned_reads.combine(downstream_reference_bundle)
        alignments = BWA_MEM2(bwa_inputs)
        sorted_bams = SORT_INDEX(alignments)
        alignment_metadata_result = ALIGNMENT_METADATA(sorted_bams.sorted_bams, reference.name, params.mt_contig)
        input_alignment_records = alignment_metadata_result[0]
        standard_alignment = input_alignment_records.map { item ->
            def fields = item[3].text.trim().split('\\n')[1].split('\t', -1)
            tuple(item[0], 'FASTQ', item[1], item[2], reference.name, fields[10], item[3])
        }
        input_metadata = input_alignment_records.map { item -> tuple(item[0], item[3]) }
    } else {
        prealigned_inputs = samples_ch
        prealigned = VALIDATE_PREALIGNED(prealigned_inputs)
        standard_alignment = prealigned.standard_alignment.map { sample, mode, bam, index, ref, metadata ->
            def fields = metadata.readLines().findAll { it.trim() }[1].split('\t', -1)
            tuple(sample, mode, bam, index, ref.toString(), fields[10], metadata)
        }
        input_metadata = prealigned.standard_alignment.map { sample, mode, bam, index, ref, metadata -> tuple(sample, metadata) }
        downstream_reference = PREPARE_DOWNSTREAM_REFERENCE(Channel.value(reference))
    }

    standardized_contract = FINALIZE_ALIGNMENT_CONTRACT(standard_alignment)
    sorted_bams = standardized_contract.contract.map { sample, mode, bam, index, ref, context, metadata -> tuple(sample, bam, index) }
    whole_flagstat = WHOLE_FLAGSTAT(sorted_bams, '')
    prealigned_flagstat = input_mode == 'FASTQ' ? Channel.empty() : INPUT_FLAGSTAT(sorted_bams, '.input')

    mt_bams = EXTRACT_MTDNA(sorted_bams, params.mt_contig)
    mtdna_flagstat = MTDNA_FLAGSTAT(mt_bams, '.mtDNA')

    mt_coverage = MTDNA_COVERAGE(mt_bams, params.mt_contig, params.coverage_mapq)
    autosomal_regions = SELECT_AUTOSOMAL_CONTIGS(sorted_bams, params.autosomal_contigs ?: '')
    autosomal_coverage = AUTOSOMAL_COVERAGE(autosomal_regions, params.coverage_mapq)
    metrics_inputs = mt_coverage.join(autosomal_coverage)
    mtdna_metrics = COMBINE_MTDNA_METRICS(metrics_inputs, params.coverage_mapq)
    reporting_depth = REPORTING_DEPTH(mt_bams, params.mt_contig, params.coverage_mapq)

    mt_reference_bundle = PREPARE_MT_REFERENCES(Channel.value(reference), params.mt_contig, params.shift_offset)
    prepared_mt_reads = PREPARE_MT_READS(mt_bams)
    mt_realign_inputs = prepared_mt_reads.combine(mt_reference_bundle)
    standard_mt_sam = REALIGN_MT_STANDARD(mt_realign_inputs, 'standard')
    shifted_mt_sam = REALIGN_MT_SHIFTED(mt_realign_inputs, 'shifted')
    standard_mt_bams = SORT_REALIGN_MT_STANDARD(standard_mt_sam)
    shifted_mt_bams = SORT_REALIGN_MT_SHIFTED(shifted_mt_sam)
    standard_mt_flagstat = REALIGN_MT_STANDARD_FLAGSTAT(standard_mt_bams)
    shifted_mt_flagstat = REALIGN_MT_SHIFTED_FLAGSTAT(shifted_mt_bams)

    gatk_reference_bundle = PREPARE_GATK_MT_REFERENCES(mt_reference_bundle)
    standard_raw = GATK_MUTECT2_STANDARD(standard_mt_bams, gatk_reference_bundle, COMBINE_MTDNA_METRICS.out, 'standard')
    shifted_raw = GATK_MUTECT2_SHIFTED(shifted_mt_bams, gatk_reference_bundle, COMBINE_MTDNA_METRICS.out, 'shifted')
    lifted_shifted = LIFTOVER_MT_VCF(shifted_raw.raw.map { sample, rep, vcf, index, stats -> tuple(sample, rep, vcf, index) }, gatk_reference_bundle)
    standard_tuple = standard_raw.raw.map { sample, rep, vcf, index, stats -> tuple(sample, vcf, index, stats) }
    lifted_with_stats = lifted_shifted.lifted
        .combine(shifted_raw.raw.map { sample, rep, vcf, index, stats -> tuple(sample, stats) })
        .map { sample, lifted_vcf, lifted_index, stats_sample, lifted_stats -> tuple(sample, lifted_vcf, lifted_index, lifted_stats) }
    canonical_inputs = standard_tuple.combine(lifted_with_stats)
        .map { sample, standard_vcf, standard_index, standard_stats, lifted_sample, lifted_vcf, lifted_index, lifted_stats -> tuple(sample, standard_vcf, standard_index, standard_stats, lifted_vcf, lifted_index, lifted_stats) }
    canonical_raw = COMBINE_GATK_VCFS(canonical_inputs, gatk_reference_bundle)
    gatk_filtered = FILTER_GATK_FINAL(canonical_raw.raw, gatk_reference_bundle, COMBINE_MTDNA_METRICS.out)

    mutserve2_native = MUTSERVE2_CALL(
        standard_mt_bams,
        gatk_reference_bundle,
        params.mt_contig,
        params.mutserve2_heteroplasmy_threshold
    )
    mutserve2_canonical = NORMALIZE_MUTSERVE2(mutserve2_native.native_calls, gatk_reference_bundle, params.mt_contig)

    gatk_interpretation_input = gatk_filtered.filtered.map { sample, vcf, index ->
        tuple(sample, 'gatk', vcf)
    }
    mutserve_interpretation_input = mutserve2_canonical.canonical.map { sample, vcf, index, calls ->
        tuple(sample, 'mutserve2', vcf)
    }

    gatk_haplogroups = HAPLOGREP_GATK(gatk_interpretation_input, params.haplogrep_tree, params.haplogrep_het_level)
    mutserve_haplogroups = HAPLOGREP_MUTSERVE2(mutserve_interpretation_input, params.haplogrep_tree, params.haplogrep_het_level)
    haplogroup_summary_inputs = gatk_haplogroups.assignment
        .map { sample, caller, standardized, native_output -> tuple(sample, standardized) }
        .join(mutserve_haplogroups.assignment.map { sample, caller, standardized, native_output -> tuple(sample, standardized) })
    combined_haplogroups = COMBINE_HAPLOGROUPS(haplogroup_summary_inputs, params.haplogrep_tree)

    gatk_contamination = HAPLOCHECK_GATK(gatk_interpretation_input)
    mutserve_contamination = HAPLOCHECK_MUTSERVE2(mutserve_interpretation_input)
    contamination_summary_inputs = gatk_contamination.assessment
        .map { sample, caller, assessment -> tuple(sample, assessment) }
        .join(mutserve_contamination.assessment.map { sample, caller, assessment -> tuple(sample, assessment) })
    combined_contamination = COMBINE_CONTAMINATION(contamination_summary_inputs)

    if (params.haplogroup_functional_fixture) {
        functional_fixture = file(params.haplogroup_functional_fixture)
        functional_expected = file(params.haplogroup_functional_expected)
        functional_input = Channel.of(tuple('NA12874', 'functional_fixture', functional_fixture))
        functional_assignment = HAPLOGREP_FUNCTIONAL(functional_input, params.haplogrep_tree, params.haplogrep_het_level)
        VALIDATE_HAPLOGREP_FUNCTIONAL(functional_assignment.assignment, functional_expected)
    }

    comparison_truth = params.fixture_truth \
        ? file(params.fixture_truth) \
        : file("${projectDir}/assets/no_fixture_truth.tsv")
    caller_inputs = canonical_raw.raw
        .join(gatk_filtered.filtered)
        .join(mutserve2_canonical.canonical)
    caller_comparison = COMPARE_CALLERS(caller_inputs, comparison_truth)

    if (params.fixture_truth) {
        fixture_truth = file(params.fixture_truth)
        VALIDATE_GATK_FIXTURE(canonical_raw.raw, fixture_truth)
        VALIDATE_MUTSERVE2_FIXTURE(mutserve2_canonical.canonical, fixture_truth)
    }

    no_numt_annotation_evidence = file("${projectDir}/assets/no_numt_variant_evidence.tsv")
    if (params.numt_strategy == 'comparison') {
        numt_capable_bams = standardized_contract.contract
            .filter { sample, mode, bam, index, ref, context, metadata -> context == 'FULL_REFERENCE' }
            .map { sample, mode, bam, index, ref, context, metadata -> tuple(sample, bam, index) }
        mapping_filtered = NUMT_MAPPING_FILTER(numt_capable_bams, params.mt_contig, params.numt_mapq)
        strategy_b_fullref = mapping_filtered.filtered_bam.map { sample, bam, index ->
            tuple("${sample}_strategy_b", bam, index)
        }
        strategy_b_reads = PREPARE_MT_READS_B(strategy_b_fullref)
        strategy_b_realign_inputs = strategy_b_reads.combine(mt_reference_bundle)
        strategy_b_standard_sam = REALIGN_MT_STANDARD_B(strategy_b_realign_inputs, 'standard')
        strategy_b_shifted_sam = REALIGN_MT_SHIFTED_B(strategy_b_realign_inputs, 'shifted')
        strategy_b_standard_bams = SORT_REALIGN_MT_STANDARD_B(strategy_b_standard_sam)
        strategy_b_shifted_bams = SORT_REALIGN_MT_SHIFTED_B(strategy_b_shifted_sam)

        strategy_b_standard_raw = GATK_MUTECT2_STANDARD_B(strategy_b_standard_bams, gatk_reference_bundle, COMBINE_MTDNA_METRICS.out, 'standard')
        strategy_b_shifted_raw = GATK_MUTECT2_SHIFTED_B(strategy_b_shifted_bams, gatk_reference_bundle, COMBINE_MTDNA_METRICS.out, 'shifted')
        strategy_b_lifted = LIFTOVER_MT_VCF_B(
            strategy_b_shifted_raw.raw.map { sample, rep, vcf, index, stats -> tuple(sample, rep, vcf, index) },
            gatk_reference_bundle
        )
        strategy_b_standard_tuple = strategy_b_standard_raw.raw.map { sample, rep, vcf, index, stats -> tuple(sample, vcf, index, stats) }
        strategy_b_lifted_with_stats = strategy_b_lifted.lifted
            .combine(strategy_b_shifted_raw.raw.map { sample, rep, vcf, index, stats -> tuple(sample, stats) })
            .map { sample, lifted_vcf, lifted_index, stats_sample, lifted_stats -> tuple(sample, lifted_vcf, lifted_index, lifted_stats) }
        strategy_b_canonical_inputs = strategy_b_standard_tuple.combine(strategy_b_lifted_with_stats)
            .map { sample, standard_vcf, standard_index, standard_stats, lifted_sample, lifted_vcf, lifted_index, lifted_stats -> tuple(sample, standard_vcf, standard_index, standard_stats, lifted_vcf, lifted_index, lifted_stats) }
        strategy_b_canonical = COMBINE_GATK_VCFS_B(strategy_b_canonical_inputs, gatk_reference_bundle)
        strategy_b_gatk_filtered = FILTER_GATK_FINAL_B(strategy_b_canonical.raw, gatk_reference_bundle, COMBINE_MTDNA_METRICS.out)

        strategy_b_mutserve_native = MUTSERVE2_CALL_B(
            strategy_b_standard_bams,
            gatk_reference_bundle,
            params.mt_contig,
            params.mutserve2_heteroplasmy_threshold
        )
        strategy_b_mutserve = NORMALIZE_MUTSERVE2_B(strategy_b_mutserve_native.native_calls, gatk_reference_bundle, params.mt_contig)

        strategy_b_gatk_for_join = strategy_b_gatk_filtered.filtered.map { sample, vcf, index ->
            tuple(sample.replaceFirst(/_strategy_b$/, ''), vcf, index)
        }
        strategy_b_mutserve_for_join = strategy_b_mutserve.canonical.map { sample, vcf, index, calls ->
            tuple(sample.replaceFirst(/_strategy_b$/, ''), vcf, index, calls)
        }

        numt_truth = params.numt_fixture_truth \
            ? file(params.numt_fixture_truth) \
            : file("${projectDir}/assets/no_numt_truth.tsv")
        numt_inputs = numt_capable_bams
            .join(standard_mt_bams)
            .join(gatk_filtered.filtered)
            .join(mutserve2_canonical.canonical)
            .join(strategy_b_gatk_for_join)
            .join(strategy_b_mutserve_for_join)
        NUMT_EVIDENCE(numt_inputs, numt_truth, params.mt_contig, params.numt_mapq)
        evaluated_numt_evidence = NUMT_EVIDENCE.out[1].map { evidence ->
            tuple(evidence.baseName.replaceFirst(/\.numt_variant_evidence$/, ''), evidence)
        }
        fallback_numt_evidence = standardized_contract.contract.map { sample, mode, bam, index, ref, context, metadata ->
            tuple(sample, no_numt_annotation_evidence)
        }
        numt_evidence_for_annotation = evaluated_numt_evidence.mix(fallback_numt_evidence)
            .groupTuple()
            .map { sample, evidence_files ->
                def evaluated = evidence_files.find { it.baseName.endsWith('.numt_variant_evidence') }
                tuple(sample, evaluated ?: evidence_files[0])
            }

        gatk_audit_inputs = canonical_raw.raw
            .join(gatk_filtered.filtered)
            .join(standard_mt_bams)
        GATK_NUMT_AUDIT(gatk_audit_inputs)

        if (params.numt_fixture_design) {
            fixture_design = file(params.numt_fixture_design)
            preflight_inputs = sorted_bams
                .join(mt_bams)
                .join(standard_mt_bams)
                .join(shifted_mt_bams)
            NUMT_FIXTURE_PREFLIGHT(preflight_inputs, fixture_design, params.mt_contig)
        }
    } else {
        numt_evidence_for_annotation = gatk_filtered.filtered.map { sample, vcf, index ->
            tuple(sample, no_numt_annotation_evidence)
        }
    }

    canonical_mt_reference = file("${projectDir}/resources/mitochondrial/NC_012920.1.fa")
    canonical_mt_features = file("${projectDir}/resources/mitochondrial/NC_012920.1.features.tsv")
    mitoarc_report_icon = file("${projectDir}/assets/images/mitoarc-icon.svg")

    gatk_annotation_inputs = gatk_filtered.filtered
        .map { sample, vcf, index -> tuple(sample, 'gatk', vcf) }
        .join(numt_evidence_for_annotation)
        .map { sample, caller, vcf, evidence -> tuple(sample, caller, vcf, evidence) }
    mutserve_annotation_inputs = mutserve2_canonical.canonical
        .map { sample, vcf, index, calls -> tuple(sample, 'mutserve2', vcf) }
        .join(numt_evidence_for_annotation)
        .map { sample, caller, vcf, evidence -> tuple(sample, caller, vcf, evidence) }

    gatk_annotations = ANNOTATE_GATK(
        gatk_annotation_inputs, gatk_reference_bundle, canonical_mt_reference,
        canonical_mt_features, params.homoplasmy_threshold
    )
    mutserve_annotations = ANNOTATE_MUTSERVE2(
        mutserve_annotation_inputs, gatk_reference_bundle, canonical_mt_reference,
        canonical_mt_features, params.homoplasmy_threshold
    )
    annotation_summary_inputs = gatk_annotations.annotation
        .map { sample, caller, annotated, summary -> tuple(sample, summary) }
        .join(mutserve_annotations.annotation.map { sample, caller, annotated, summary -> tuple(sample, summary) })
    combined_annotation_summary = COMBINE_ANNOTATION_SUMMARIES(annotation_summary_inputs)

    gatk_consensus = CONSENSUS_GATK(
        gatk_filtered.filtered.map { sample, vcf, index -> tuple(sample, 'gatk', vcf) },
        gatk_reference_bundle, canonical_mt_reference, params.homoplasmy_threshold,
        params.consensus_heteroplasmy_mode
    )
    mutserve2_consensus = CONSENSUS_MUTSERVE2(
        mutserve2_canonical.canonical.map { sample, vcf, index, calls -> tuple(sample, 'mutserve2', vcf) },
        gatk_reference_bundle, canonical_mt_reference, params.homoplasmy_threshold,
        params.consensus_heteroplasmy_mode
    )

    if (params.annotation_functional_fixture) {
        VALIDATE_ANNOTATION_FUNCTIONAL(
            file(params.annotation_functional_fixture), file(params.annotation_functional_expected),
            canonical_mt_reference, canonical_mt_features, no_numt_annotation_evidence,
            params.homoplasmy_threshold
        )
        REPORTING_FUNCTIONAL_FIXTURE(
            file(params.annotation_functional_fixture), file(params.annotation_functional_expected),
            canonical_mt_reference, canonical_mt_features, no_numt_annotation_evidence,
            mitoarc_report_icon,
            params.coverage_mapq, params.homoplasmy_threshold,
            params.mutserve2_heteroplasmy_threshold, workflow.manifest.version ?: 'NA',
            workflow.nextflow.version.toString()
        )
    }

    metrics_for_reporting = mtdna_metrics.map { metric ->
        tuple(metric.baseName.replaceFirst(/\.mtdna_metrics$/, ''), metric)
    }
    caller_comparison_for_reporting = caller_comparison[0]
        .map { comparison -> tuple(comparison.baseName.replaceFirst(/\.caller_comparison$/, ''), comparison) }
        .join(caller_comparison[1].map { summary ->
            tuple(summary.baseName.replaceFirst(/\.caller_summary$/, ''), summary)
        })
    haplogroups_for_reporting = combined_haplogroups[0]
        .map { haplogroups -> tuple(haplogroups.baseName.replaceFirst(/\.haplogroups$/, ''), haplogroups) }
        .join(combined_haplogroups[1].map { comparison ->
            tuple(comparison.baseName.replaceFirst(/\.haplogroup_comparison$/, ''), comparison)
        })
    contamination_for_reporting = combined_contamination.map { contamination ->
        tuple(contamination.baseName.replaceFirst(/\.contamination$/, ''), contamination)
    }
    annotation_summary_for_reporting = combined_annotation_summary.map { summary ->
        tuple(summary.baseName.replaceFirst(/\.annotation_summary$/, ''), summary)
    }
    annotations_for_reporting = gatk_annotations.annotation
        .map { sample, caller, annotated, summary -> tuple(sample, annotated) }
        .join(mutserve_annotations.annotation.map { sample, caller, annotated, summary -> tuple(sample, annotated) })
    consensus_for_reporting = gatk_consensus.consensus
        .map { sample, caller, fasta, validation -> tuple(sample, validation) }
        .join(mutserve2_consensus.consensus.map { sample, caller, fasta, validation -> tuple(sample, validation) })

    reporting_inputs = reporting_depth.depth
        .join(metrics_for_reporting)
        .join(caller_comparison_for_reporting)
        .join(haplogroups_for_reporting)
        .join(contamination_for_reporting)
        .join(annotation_summary_for_reporting)
        .join(annotations_for_reporting)
        .join(consensus_for_reporting)
        .join(input_metadata.map { sample, metadata -> tuple(sample, metadata) })
        .map { sample, depth, metrics, source_comparison, caller_summary,
              haplogroups, haplogroup_comparison, contamination, annotation_summary,
              gatk_annotation, mutserve2_annotation, gatk_consensus_validation,
              mutserve2_consensus_validation, input_validation ->
            tuple(sample, depth, metrics, source_comparison, caller_summary, haplogroups,
                  haplogroup_comparison, contamination, annotation_summary, gatk_annotation,
                  mutserve2_annotation, gatk_consensus_validation,
                  mutserve2_consensus_validation, input_validation)
        }

    reporting_contract = PREPARE_REPORTING_DATA(
        reporting_inputs, canonical_mt_features, reference.name, params.mt_contig,
        params.coverage_mapq, params.homoplasmy_threshold,
        params.mutserve2_heteroplasmy_threshold, params.numt_strategy,
        workflow.manifest.version ?: 'NA', workflow.nextflow.version.toString()
    )
    reporting_plots = GENERATE_REPORT_PLOTS(reporting_contract.contract)
    report_render_inputs = reporting_contract.contract
        .map { sample, summary, summary_json, plot_data, features, variants, comparison ->
            tuple(sample, summary, variants, comparison)
        }
        .join(reporting_plots.plots.map { sample, svgs, pngs, manifest -> tuple(sample, svgs, pngs) })
    integrated_report = RENDER_MTDNA_REPORT(report_render_inputs, mitoarc_report_icon)

    reporting_validation_sources = reporting_inputs.map {
            sample, depth, metrics, source_comparison, caller_summary, haplogroups,
        haplogroup_comparison, contamination, annotation_summary, gatk_annotation,
        mutserve2_annotation, gatk_consensus_validation, mutserve2_consensus_validation, input_validation ->
            tuple(sample, metrics, annotation_summary, gatk_annotation, mutserve2_annotation, source_comparison)
    }
    reporting_validation_inputs = reporting_contract.contract
        .map { sample, summary, summary_json, plot_data, features, variants, comparison ->
            tuple(sample, summary, plot_data, features, variants, comparison)
        }
        .join(reporting_validation_sources)
        .join(reporting_plots.plots)
        .join(integrated_report.report)
    VALIDATE_REPORTING(reporting_validation_inputs)

    if (input_mode == 'FASTQ') {
        fastqc_multiqc_files = fastqc_reports.flatMap { sample, reads1, reads2, html_reports, zip_reports ->
            def reports = zip_reports instanceof Collection ? zip_reports : [zip_reports]
            reports
        }
        multiqc_files = fastqc_multiqc_files.mix(FASTP.out[2])
    } else {
        multiqc_files = Channel.empty()
    }
    multiqc_files = multiqc_files
        .mix(whole_flagstat, prealigned_flagstat, mtdna_flagstat, standard_mt_flagstat, shifted_mt_flagstat)
        .collect()
    multiqc_result = MULTIQC_REPORT(multiqc_files)
    VALIDATE_MULTIQC_REPORT(multiqc_result[0], multiqc_result[1])
}
