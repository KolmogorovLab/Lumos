#!/usr/bin/env nextflow
nextflow.enable.dsl     = 2
//nextflow.preview.output = true   //uncomment for for 25.04.2


include { alignMinimap2; callClair3; phaseLongphase; deepsomaticTumorOnly;
          modkitDMR; modkitPileupAllele; modkitPileup;
          severusTumorOnly; haplotagWhatshap; wakhanCNA; wakhanHapcorrect } from "./processes/processes.nf"
include { modkitStats as modkitStats  } from "./processes/processes.nf"
include { modkitStats as modkitStats2 } from "./processes/processes.nf"
include { modkitStats as modkitStats3 } from "./processes/processes.nf"
include { fastaIndex as fastaIndex } from "./processes/processes.nf"


// ---------- default user params ----------
params.reads_tumor = null
params.reference = null
params.vntr = null
params.sv_pon = null
params.clair3_model = null
params.cpgs = null
params.cosmic = null

params.mode = 'all'         // all | sv_cna | sv_cna_dmr

params.aligned_input = 'false'        // 'true' | 'false'
params.aligned_tumor = null
params.aligned_tumor_bai = null

// ---------- subworkflow with alignment toggle + modes ----------
workflow tumorOnlyOntWorkflow {

  take:
    reads_tumor               // (sample, fastq) tuples (used only if alignment == 'true')
    reference           // fasta
    alignedTumor       // BAM (used only if alignment == 'false')
    alignedTumorBai       // BAI (used only if alignment == 'false')
    vntrAnnotation
    svPanelOfNormals
    clair3Model
    cpgs
    cosmic

  main:
    def RUN_ALL     = (params.mode == 'all')
    def RUN_SV_CNA  = (params.mode in ['all','sv_cna','sv_cna_dmr'])
    def RUN_DMR     = (params.mode in ['all','sv_cna_dmr'])
    def RUN_DEEPSOM = (params.mode == 'all')

    // Unify BAM/BAI/REF_IDX depending on alignment mode
    def alignTrue = params.aligned_input.toString().toLowerCase() == 'false'

    // Defaults
    def bamCh
    def baiCh
    def refIdxCh

    if (alignTrue) {
      // Align from reads
      alignMinimap2(reference, reads_tumor.collect())
      bamCh    = alignMinimap2.out.bam
      baiCh    = alignMinimap2.out.bam_idx
      refIdxCh = alignMinimap2.out.ref_idx
    }
    else {
      // Use pre-aligned BAM/BAI; ensure we have a fasta index
      fastaIndex(reference)
      bamCh    = alignedTumor
      baiCh    = alignedTumorBai
      refIdxCh = fastaIndex.out.ref_idx
    }

    // 2) Clair3 (requires BAM/BAI/REF/REF_IDX)
    callClair3(
      bamCh,
      baiCh,
      reference,
      refIdxCh,
      clair3Model
    )

    // 3) Longphase (requires BAM/BAI/REF/REF_IDX + VCF)
    phaseLongphase(
      bamCh,
      baiCh,
      reference,
      refIdxCh,
      callClair3.out.vcf
    )

    // 4) Wakhan hap-correct
    wakhanHapcorrect(
      bamCh,
      baiCh,
      reference,
      phaseLongphase.out.phasedVcf
    )

    // 5) Whatshap haplotag
    haplotagWhatshap(
      reference,
      refIdxCh,
      wakhanHapcorrect.out.rephasedVcf,
      bamCh,
      baiCh
    )

    // Prepare default (maybe-empty) channels for conditional steps
    def severusSomaticVcfCh = Channel.empty()
    def severusFullOutputCh = Channel.empty()
    def wakhanDataOutputCh = Channel.empty()
    def wakhanFullOutputCh  = Channel.empty()
    def hp1Ch = Channel.empty(); def hp2Ch = Channel.empty()
    def pileCh = Channel.empty(); def pile2Ch = Channel.empty(); def dmrCh = Channel.empty()
    def s1Ch = Channel.empty();  def s2Ch  = Channel.empty(); def s3Ch = Channel.empty()
    def deepSomCh = Channel.empty()

    // 6) Severus + 7) WakhanCNA (if SV/CNA enabled)
    if (RUN_SV_CNA) {
      severusTumorOnly(
        haplotagWhatshap.out.bam,
        haplotagWhatshap.out.bam_idx,
        wakhanHapcorrect.out.rephasedVcf,
        vntrAnnotation,
        svPanelOfNormals
      )

      // keep the 6th input from wakhanHapcorrect
      wakhanCNA(
        bamCh,
        baiCh,
        reference,
        wakhanHapcorrect.out.rephasedVcf,
        severusTumorOnly.out.severusSomaticVcf,
        wakhanHapcorrect.out.wakhanHPOutput,
        cosmic
      )

      severusSomaticVcfCh = severusTumorOnly.out.severusSomaticVcf
      severusFullOutputCh = severusTumorOnly.out.severusFullOutput
      wakhanFullOutputCh  = wakhanCNA.out.wakhanOutput
      wakhanDataOutputCh = wakhanHapcorrect.out.wakhanHPOutput
    }

    // 8) Modkit (DMR/Stats) if enabled
    if (RUN_DMR) {
      modkitPileupAllele(
        haplotagWhatshap.out.bam,
        haplotagWhatshap.out.bam_idx,
        reference,
        refIdxCh
      )

      modkitPileup(
        haplotagWhatshap.out.bam,
        haplotagWhatshap.out.bam_idx,
        reference,
        refIdxCh,
	cpgs
      )

      modkitDMR(
        modkitPileupAllele.out.HP1bed,
        modkitPileupAllele.out.HP1bed_idx,
        modkitPileupAllele.out.HP2bed,
        modkitPileupAllele.out.HP2bed_idx,
        reference,
        refIdxCh,
        cpgs
      )

      modkitStats(
        modkitPileupAllele.out.HP1bed,
        modkitPileupAllele.out.HP1bed_idx,
        reference,
        refIdxCh,
        cpgs
      )

      modkitStats2(
        modkitPileupAllele.out.HP2bed,
        modkitPileupAllele.out.HP2bed_idx,
        reference,
        refIdxCh,
        cpgs
      )

      modkitStats3(
        modkitPileup.out.pileupbed,
        modkitPileup.out.pileupbed_idx,
        reference,
        refIdxCh,
        cpgs
      )

      hp1Ch   = modkitPileupAllele.out.HP1bed
      hp2Ch   = modkitPileupAllele.out.HP2bed
      pileCh  = modkitPileup.out.pileupbed
      pile2Ch = modkitPileup.out.pileupbed_subset
      dmrCh   = modkitDMR.out.DMRbed
      s1Ch    = modkitStats.out.stats
      s2Ch    = modkitStats2.out.stats
      s3Ch    = modkitStats3.out.stats
    }

    // 9) DeepSomatic only in 'all'
    if (RUN_DEEPSOM) {
      deepsomaticTumorOnly(
        bamCh,
        baiCh,
        reference,
        refIdxCh
      )
      deepSomCh = deepsomaticTumorOnly.out.deepsomaticOutput
    }

  emit:
    // phasing / hap-correction / haplotag
    phasedVcf              = phaseLongphase.out.phasedVcf
    rephasedVcf            = wakhanHapcorrect.out.rephasedVcf
    haplotaggedBam         = haplotagWhatshap.out.bam
    haplotaggedBamidx      = haplotagWhatshap.out.bam_idx

    // SV/CNA
    severusFullOutput      = severusFullOutputCh
    wakhanFullOutput       = wakhanFullOutputCh
    wakhanDataOutput       = wakhanDataOutputCh

    // DeepSomatic
    deepsomaticOutput      = deepSomCh

    // Methylation
    modkitPileupAlleleBED1 = hp1Ch
    modkitPileupAlleleBED2 = hp2Ch
    modkitPileupOut        = pileCh
    modkitPileup2Out       = pile2Ch
    modkitDMROut           = dmrCh
    modkitStatsOut         = s1Ch
    modkitStats2Out        = s2Ch
    modkitStats3Out        = s3Ch
}

// ---------- entry workflow (publish/output) ----------
workflow {

  main:
    // Arg checks depend on alignment mode
    def missing = []
    if (!params.reference)    missing << '--reference'
    if (!params.vntr)         missing << '--vntr'
    if (!params.sv_pon)       missing << '--sv_pon'
    if (!params.clair3_model) missing << '--clair3_model'
    if (!params.cpgs)         missing << '--cpgs'

    def alignTrue = params.aligned_input.toString().toLowerCase() == 'false'
    if (alignTrue) {
      if (!params.reads_tumor) missing << '--reads_tumor'
    } else {
      if (!params.aligned_tumor) missing << '--aligned_tumor'
      if (!params.aligned_tumor_bai) missing << '--aligned_tumor_bai'
    }
    if (missing) error "Missing required arguments: ${missing.join(', ')}"

    log.info "Mode: ${params.mode} | Pre-aligned input: ${params.aligned_input}"

    
    // Build channels for whichever path we use
    reads_ch = alignTrue
       ? Channel.fromPath(params.reads_tumor.split(" ").toList(), checkIfExists: true)
       : Channel.empty()
    reads_ch.view{it -> "Input reads: $it"}
    
    pre_bam_ch = !alignTrue ? Channel.fromPath(params.aligned_tumor, checkIfExists:true) : Channel.empty()
    pre_bai_ch = !alignTrue ? Channel.fromPath(params.aligned_tumor_bai, checkIfExists:true) : Channel.empty()

    pre_bam_ch.view { "Tumor BAM: $it" }
    pre_bai_ch.view { "Tumor BAI: $it" }

    ref_ch    = Channel.fromPath(params.reference,    checkIfExists:true)
    vntr_ch   = Channel.fromPath(params.vntr,         checkIfExists:true)
    svpon_ch  = Channel.fromPath(params.sv_pon,       checkIfExists:true)
    clair3_ch = Channel.fromPath(params.clair3_model, checkIfExists:true)
    cpgs_ch   = Channel.fromPath(params.cpgs,         checkIfExists:true)

    cosmic_ch = params.cosmic 
        ? Channel.fromPath(params.cosmic, checkIfExists:true)
        : Channel.fromPath(file('NOFILE'))

    out = tumorOnlyOntWorkflow(
      reads_ch, ref_ch, pre_bam_ch, pre_bai_ch, vntr_ch, svpon_ch, clair3_ch, cpgs_ch, cosmic_ch
    )

    publish:
        phasedVcf              = out.phasedVcf
        rephasedVcf            = out.rephasedVcf
        haplotaggedBam         = out.haplotaggedBam
        haplotaggedBamidx      = out.haplotaggedBamidx
        severusFullOutput      = out.severusFullOutput
        wakhanFullOutput       = out.wakhanFullOutput
        wakhanDataOutput       = out.wakhanDataOutput
        deepsomaticOutput      = out.deepsomaticOutput
        modkitPileupAlleleBED1 = out.modkitPileupAlleleBED1
        modkitPileupAlleleBED2 = out.modkitPileupAlleleBED2
        modkitPileupOut        = out.modkitPileupOut
        modkitPileupClean      = out.modkitPileup2Out
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
    wakhanDataOutput       { path 'wakhan'  }
    deepsomaticOutput      { path 'deepsomatic' }
    modkitPileupAlleleBED1 { path 'methylation' }
    modkitPileupAlleleBED2 { path 'methylation' }
    modkitPileupOut        { path 'methylation' }
    modkitPileupClean      { path 'methylation' }
    modkitDMROut           { path 'methylation' }
    modkitStatsOut         { path 'methylation' }
    modkitStats2Out        { path 'methylation' }
    modkitStats3Out        { path 'methylation' }
  }


