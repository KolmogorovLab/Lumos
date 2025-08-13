#!/usr/bin/env nextflow
nextflow.enable.dsl     = 2
nextflow.preview.output = true   // required for entry-workflow publish/output

// ---------- modules (assumes these exist with matching IO) ----------
include { alignMinimap2; callClair3; phaseLongphase; deepsomaticTumorOnly;
          modkitDMR; modkitPileupAllele; modkitPileup;
          severusTumorOnly; haplotagWhatshap; wakhanCNA; wakhanHapcorrect } from "./processes/processes.nf"
include { modkitStats as modkitStats  } from "./processes/processes.nf"
include { modkitStats as modkitStats2 } from "./processes/processes.nf"
include { modkitStats as modkitStats3 } from "./processes/processes.nf"

// ---------- user params ----------
params.reads        = params.reads        ?: null
params.reference    = params.reference    ?: null
params.vntr         = params.vntr         ?: null
params.sv_pon       = params.sv_pon       ?: null
params.clair3_model = params.clair3_model ?: null
params.cpgs         = params.cpgs         ?: null
params.mode         = params.mode         ?: 'all'   // all | sv_cna | sv_cna_dmr
params.alignment    = params.alignment    ?: 'true'  // (placeholder)

// ---------- subworkflow with conditional steps ----------
workflow tumorOnlyOntWorkflow {

  take:
    reads
    reference
    vntrAnnotation
    svPanelOfNormals
    clair3Model
    cpgs

  main:
    // mode flags
    def RUN_ALL        = (params.mode == 'all')
    def RUN_SV_CNA     = (params.mode in ['all','sv_cna','sv_cna_dmr'])
    def RUN_DMR        = (params.mode in ['all','sv_cna_dmr'])
    def RUN_DEEPSOM    = (params.mode == 'all')

    // 1) Align (always, as written)
    alignMinimap2(reference, reads.collect())

    // 2) Clair3 (always; longphase needs VCF)
    callClair3(
      alignMinimap2.out.bam,
      alignMinimap2.out.bam_idx,
      reference,
      alignMinimap2.out.ref_idx,
      clair3Model
    )

    // 3) Longphase (always)
    phaseLongphase(
      alignMinimap2.out.bam,
      alignMinimap2.out.bam_idx,
      reference,
      alignMinimap2.out.ref_idx,
      callClair3.out.vcf
    )

    // 4) Wakhan hap-correct (always)
    wakhanHapcorrect(
      alignMinimap2.out.bam,
      alignMinimap2.out.bam_idx,
      reference,
      phaseLongphase.out.phasedVcf
    )

    // 5) Whatshap haplotag (always)
    haplotagWhatshap(
      reference,
      alignMinimap2.out.ref_idx,
      wakhanHapcorrect.out.rephasedVcf,
      alignMinimap2.out.bam,
      alignMinimap2.out.bam_idx
    )

    // Defaults for emits in case we skip things
    def severusSomaticVcf    = Channel.empty()
    def severusFullOutputCh  = Channel.empty()
    def wakhanFullOutputCh   = Channel.empty()
    def modkitHP1Ch          = Channel.empty()
    def modkitHP2Ch          = Channel.empty()
    def modkitPileCh         = Channel.empty()
    def dmrCh                = Channel.empty()
    def stats1Ch             = Channel.empty()
    def stats2Ch             = Channel.empty()
    def stats3Ch             = Channel.empty()
    def deepSomCh            = Channel.empty()

    // 6) Severus + 7) WakhanCNA (conditional on SV/CNA modes)
    if ( RUN_SV_CNA ) {
      severusTumorOnly(
        haplotagWhatshap.out.bam,
        haplotagWhatshap.out.bam_idx,
        wakhanHapcorrect.out.rephasedVcf,
        vntrAnnotation,
        svPanelOfNormals
      )

      // keep 6th input from wakhanHapcorrect
      wakhanCNA(
        alignMinimap2.out.bam,
        alignMinimap2.out.bam_idx,
        reference,
        wakhanHapcorrect.out.rephasedVcf,
        severusTumorOnly.out.severusSomaticVcf,
        wakhanHapcorrect.out.wakhanOutput
      )

      severusSomaticVcf   = severusTumorOnly.out.severusSomaticVcf
      severusFullOutputCh = severusTumorOnly.out.severusFullOutput
      wakhanFullOutputCh  = wakhanCNA.out.wakhanOutput
    }

    // 8) Modkit (only in all or sv_cna_dmr)
    if ( RUN_DMR ) {
      modkitPileupAllele(
        haplotagWhatshap.out.bam,
        haplotagWhatshap.out.bam_idx,
        reference,
        alignMinimap2.out.ref_idx
      )

      modkitPileup(
        haplotagWhatshap.out.bam,
        haplotagWhatshap.out.bam_idx,
        reference,
        alignMinimap2.out.ref_idx
      )

      modkitDMR(
        modkitPileupAllele.out.HP1bed,
        modkitPileupAllele.out.HP1bed_idx,
        modkitPileupAllele.out.HP2bed,
        modkitPileupAllele.out.HP2bed_idx,
        reference,
        alignMinimap2.out.ref_idx,
        cpgs
      )

      modkitStats(
        modkitPileupAllele.out.HP1bed,
        modkitPileupAllele.out.HP1bed_idx,
        reference,
        alignMinimap2.out.ref_idx,
        cpgs
      )

      modkitStats2(
        modkitPileupAllele.out.HP2bed,
        modkitPileupAllele.out.HP2bed_idx,
        reference,
        alignMinimap2.out.ref_idx,
        cpgs
      )

      modkitStats3(
        modkitPileup.out.pileupbed,
        modkitPileup.out.pileupbed_idx,
        reference,
        alignMinimap2.out.ref_idx,
        cpgs
      )

      modkitHP1Ch  = modkitPileupAllele.out.HP1bed
      modkitHP2Ch  = modkitPileupAllele.out.HP2bed
      modkitPileCh = modkitPileup.out.pileupbed
      dmrCh        = modkitDMR.out.DMRbed
      stats1Ch     = modkitStats.out.stats
      stats2Ch     = modkitStats2.out.stats
      stats3Ch     = modkitStats3.out.stats
    }

    // 9) DeepSomatic (only in `all`)
    if ( RUN_DEEPSOM ) {
      deepsomaticTumorOnly(
        alignMinimap2.out.bam,
        alignMinimap2.out.bam_idx,
        reference,
        alignMinimap2.out.ref_idx
      )
      deepSomCh = deepsomaticTumorOnly.out.deepsomaticOutput
    }

  emit:
    // phasing & hap-correction / haplotag
    phasedVcf              = phaseLongphase.out.phasedVcf
    rephasedVcf            = wakhanHapcorrect.out.rephasedVcf
    haplotaggedBam         = haplotagWhatshap.out.bam
    haplotaggedBamidx      = haplotagWhatshap.out.bam_idx
    // SV/CNA (maybe empty if mode excludes)
    severusFullOutput      = severusFullOutputCh
    wakhanFullOutput       = wakhanFullOutputCh
    // DeepSomatic (maybe empty)
    deepsomaticOutput      = deepSomCh
    // Methylation (maybe empty)
    modkitPileupAlleleBED1 = modkitHP1Ch
    modkitPileupAlleleBED2 = modkitHP2Ch
    modkitPileupOut        = modkitPileCh
    modkitDMROut           = dmrCh
    modkitStatsOut         = stats1Ch
    modkitStats2Out        = stats2Ch
    modkitStats3Out        = stats3Ch
}

// ---------- entry workflow (publish/output) ----------
workflow {

  main:
    // arg checks
    def missing = []
    if (!params.reads)        missing << '--reads'
    if (!params.reference)    missing << '--reference'
    if (!params.vntr)         missing << '--vntr'
    if (!params.sv_pon)       missing << '--sv_pon'
    if (!params.clair3_model) missing << '--clair3_model'
    if (!params.cpgs)         missing << '--cpgs'
    if (missing) error "Missing required arguments: ${missing.join(', ')}"

    log.info "Mode: ${params.mode} | Alignment: ${params.alignment}"

    reads_ch = Channel.fromPath(params.reads.split(" ").toList(), checkIfExists: true)
    reads_ch.view{it -> "Input reads: $it"}
    ref_ch    = Channel.fromPath(params.reference,    checkIfExists:true)
    vntr_ch   = Channel.fromPath(params.vntr,         checkIfExists:true)
    svpon_ch  = Channel.fromPath(params.sv_pon,       checkIfExists:true)
    clair3_ch = Channel.fromPath(params.clair3_model, checkIfExists:true)
    cpgs_ch   = Channel.fromPath(params.cpgs,         checkIfExists:true)

    out = tumorOnlyOntWorkflow(reads_ch, ref_ch, vntr_ch, svpon_ch, clair3_ch, cpgs_ch)

 publish:
    // assign to workflow outputs
    phasedVcf              = out.phasedVcf
    rephasedVcf            = out.rephasedVcf
    haplotaggedBam         = out.haplotaggedBam
    haplotaggedBamidx      = out.haplotaggedBamidx
    severusFullOutput      = out.severusFullOutput
    wakhanFullOutput       = out.wakhanFullOutput
    deepsomaticOutput      = out.deepsomaticOutput
    modkitPileupAlleleBED1 = out.modkitPileupAlleleBED1
    modkitPileupAlleleBED2 = out.modkitPileupAlleleBED2
    modkitPileupOut        = out.modkitPileupOut
    modkitDMROut           = out.modkitDMROut
    modkitStatsOut         = out.modkitStatsOut
    modkitStats2Out        = out.modkitStats2Out
    modkitStats3Out        = out.modkitStats3Out
}
  output {
    phasedVcf              { path 'phased_vcf'   }
    rephasedVcf            { path 'rephased_vcf' }
    haplotaggedBam         { path 'haplotagged_bam' }
    haplotaggedBamidx      { path 'haplotagged_bam' }
    severusFullOutput      { path 'severus' }
    wakhanFullOutput       { path 'wakhan'  }
    deepsomaticOutput      { path 'deepsomatic' }
    modkitPileupAlleleBED1 { path 'methylation' }
    modkitPileupAlleleBED2 { path 'methylation' }
    modkitPileupOut        { path 'methylation' }
    modkitDMROut           { path 'methylation' }
    modkitStatsOut         { path 'methylation' }
    modkitStats2Out        { path 'methylation' }
    modkitStats3Out        { path 'methylation' }
  }

