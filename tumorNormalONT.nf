#!/usr/bin/env nextflow
nextflow.enable.dsl     = 2
nextflow.preview.output = true   


include { callClair3; phaseLongphase; severusTumorNormal; wakhanCNATN; wakhanHapcorrectTN; modkitDMR; modkitPileup; modkitPileupAllele; deepsomaticTumorNormal } from "./processes/processes.nf"
include { alignMinimap2 as alignTumor } from "./processes/processes.nf"
include { alignMinimap2 as alignNormal } from "./processes/processes.nf"
include { haplotagWhatshap as haplotagNormal } from "./processes/processes.nf"
include { haplotagWhatshap as haplotagTumor } from "./processes/processes.nf"
include { modkitStats as modkitStats  } from "./processes/processes.nf"
include { modkitStats as modkitStats2 } from "./processes/processes.nf"
include { modkitStats as modkitStats3 } from "./processes/processes.nf"


process REF_INDEX {
  tag "${reference?.name}"
  input:
    path reference
  output:
    path "${reference}.fai", emit: ref_idx
  script:
  """
  samtools faidx ${reference}
  """
}

// ---------- user params ----------
params.readsT         = params.readsT          ?: null
params.readsN         = params.readsN          ?: null   
params.prealignedBam  = params.prealignedBam   ?: null
params.prealignedBamN = params.prealignedBamN  ?: null     
params.reference      = params.reference       ?: null
params.vntr           = params.vntr            ?: null
params.clair3_model   = params.clair3_model    ?: null
params.cpgs           = params.cpgs            ?: null
params.bam            = params.bam             ?: null          // used when --alignment 'false'
params.bai            = params.bai             ?: null          // used when --alignment 'false'
params.mode           = params.mode            ?: 'all'         // all | sv_cna | sv_cna_dmr
params.alignment      = params.alignment       ?: 'true'        // 'true' | 'false'

// ---------- subworkflow with alignment toggle + modes ----------
workflow tumorNormalOntWorkflow {

  take:
    readsT               // (Tumor sample, fastq) tuples (used only if alignment == 'true')
    readsN               // (Normal sample, fastq) tuples (used only if alignment == 'true')
    prealignedBam       // BAM (used only if alignment == 'false')
    prealignedBai       // BAI (used only if alignment == 'false')
    NprealignedBam       // BAM (used only if alignment == 'false')
    NprealignedBai       // BAI (used only if alignment == 'false')
    reference           // fasta
    vntrAnnotation
    clair3Model
    cpgs

  main:
    def RUN_ALL     = (params.mode == 'all')
    def RUN_SV_CNA  = (params.mode in ['all','sv_cna','sv_cna_dmr'])
    def RUN_DMR     = (params.mode in ['all','sv_cna_dmr'])
    def RUN_DEEPSOM = (params.mode == 'all')

    // Unify BAM/BAI/REF_IDX depending on alignment mode
    Channel
      .value(params.alignment.toString().toLowerCase() == 'true')
      .set { ALIGN_FLAG }

    // Defaults
    def bamCh
    def baiCh
    def refIdxCh

    if (ALIGN_FLAG.first()) {
      // Align from reads
      alignTumor(reference, readsT.collect())
      bamCh    = alignTumor.out.bam
      baiCh    = alignTumor.out.bam_idx
      refIdxCh = alignTumor.out.ref_idx
	  
      alignNormal(reference, readsN.collect())
      NbamCh    = alignNormal.out.bam
      NbaiCh    = alignNormal.out.bam_idx
    }
    else {
      // Use pre-aligned BAM/BAI; ensure we have a fasta index
      REF_INDEX(reference)
      bamCh    = prealignedBam
      baiCh    = prealignedBai
      NbamCh    = NprealignedBam
      NbaiCh    = NprealignedBai
      refIdxCh = REF_INDEX.out.ref_idx
    }

    // 2) Clair3 (requires BAM/BAI/REF/REF_IDX)
    callClair3(
      NbamCh,
      NbaiCh,
      reference,
      refIdxCh,
      clair3Model
    )

    // 3) Longphase (requires BAM/BAI/REF/REF_IDX + VCF)
    phaseLongphase(
      NbamCh,
      NbaiCh,
      reference,
      refIdxCh,
      callClair3.out.vcf
    )

    // 4) Wakhan hap-correct
    wakhanHapcorrectTN(
      bamCh,
      baiCh,
      reference,
      phaseLongphase.out.phasedVcf
    )

    // 5) Whatshap haplotag
    haplotagNormal(
      reference,
      refIdxCh,
      wakhanHapcorrectTN.out.rephasedVcf,
      NbamCh,
      NbaiCh
    )
    haplotagTumor(
      reference,
      refIdxCh,
      wakhanHapcorrectTN.out.rephasedVcf,
      bamCh,
      baiCh
    )

    // Prepare default (maybe-empty) channels for conditional steps
    def severusSomaticVcfCh = Channel.empty()
    def severusFullOutputCh = Channel.empty()
    def wakhanDataOutputCh = Channel.empty()
    def wakhanFullOutputCh  = Channel.empty()
    def hp1Ch = Channel.empty(); def hp2Ch = Channel.empty()
    def pileCh = Channel.empty(); def dmrCh = Channel.empty()
    def s1Ch = Channel.empty();  def s2Ch  = Channel.empty(); def s3Ch = Channel.empty()
    def deepSomCh = Channel.empty()

    // 6) Severus + 7) WakhanCNA (if SV/CNA enabled)
    if (RUN_SV_CNA) {
      severusTumorNormal(
        haplotagTumor.out.bam,
        haplotagTumor.out.bam_idx,
        haplotagNormal.out.bam,
        haplotagNormal.out.bam_idx,
        wakhanHapcorrectTN.out.rephasedVcf,
        vntrAnnotation
      )

      // keep the 6th input from wakhanHapcorrect
      wakhanCNATN(
        bamCh,
        baiCh,
        reference,
        wakhanHapcorrectTN.out.rephasedVcf,
        severusTumorNormal.out.severusSomaticVcf,
        wakhanHapcorrectTN.out.wakhanHPOutput
      )

      severusSomaticVcfCh = severusTumorNormal.out.severusSomaticVcf
      severusFullOutputCh = severusTumorNormal.out.severusFullOutput
      wakhanFullOutputCh  = wakhanCNATN.out.wakhanOutput
      wakhanDataOutputCh = wakhanHapcorrectTN.out.wakhanHPOutput
    }

    // 8) Modkit (DMR/Stats) if enabled
    if (RUN_DMR) {
      modkitPileupAllele(
        haplotagTumor.out.bam,
        haplotagTumor.out.bam_idx,
        reference,
        refIdxCh
      )

      modkitPileup(
        haplotagTumor.out.bam,
        haplotagTumor.out.bam_idx,
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
    rephasedVcf            = wakhanHapcorrectTN.out.rephasedVcf
    haplotaggedBam         = haplotagTumor.out.bam
    haplotaggedBamidx      = haplotagTumor.out.bam_idx
    haplotaggedBamN        = haplotagNormal.out.bam
    haplotaggedBamidxN     = haplotagNormal.out.bam_idx

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
    if (!params.clair3_model) missing << '--clair3_model'
    if (!params.cpgs)         missing << '--cpgs'

    def alignTrue = params.alignment.toString().toLowerCase() == 'true'
    if (alignTrue) {
      if (!params.tumor_reads) missing << '--tumor_reads'
	  if (!params.normal_reads) missing << '--normal_reads'
	  
    } else {
      if (!params.tumor_bam) missing << '--tumor_bam'
      if (!params.tumor_bai) missing << '--tumor_bai'

      if (!params.normal_bam) missing << '--tumor_bam'
      if (!params.normal_bai) missing << '--tumor_bai'
    }
    if (missing) error "Missing required arguments: ${missing.join(', ')}"

    log.info "Mode: ${params.mode} | Alignment: ${params.alignment}"

    
    // Build channels for whichever path we use
    readsT_ch = alignTrue
       ? Channel.fromPath(params.tumor_reads.split(" ").toList(), checkIfExists: true)
       : Channel.empty()
    readsT_ch.view{it -> "Input reads: $it"}
	
    readsN_ch = alignTrue
       ? Channel.fromPath(params.normal_reads.split(" ").toList(), checkIfExists: true)
       : Channel.empty()
    
    pre_bam_ch = !alignTrue ? Channel.fromPath(params.bam, checkIfExists:true) : Channel.empty()
    pre_bai_ch = !alignTrue ? Channel.fromPath(params.bai, checkIfExists:true) : Channel.empty()
	
    pre_bamN_ch = !alignTrue ? Channel.fromPath(params.normal_bam, checkIfExists:true) : Channel.empty()
    pre_baiN_ch = !alignTrue ? Channel.fromPath(params.normal_bai, checkIfExists:true) : Channel.empty()

    ref_ch    = Channel.fromPath(params.reference,    checkIfExists:true)
    vntr_ch   = Channel.fromPath(params.vntr,         checkIfExists:true)
    clair3_ch = Channel.fromPath(params.clair3_model, checkIfExists:true)
    cpgs_ch   = Channel.fromPath(params.cpgs,         checkIfExists:true)

    out = tumorNormalOntWorkflow(
      readsT_ch, readsN_ch, pre_bam_ch, pre_bai_ch, pre_bamN_ch, pre_baiN_ch, ref_ch, vntr_ch, clair3_ch, cpgs_ch
    )

    publish:
	phasedVcf              = out.phasedVcf
	rephasedVcf            = out.rephasedVcf
	haplotaggedBam         = out.haplotaggedBam
	haplotaggedBamidx      = out.haplotaggedBamidx
	haplotaggedBamN        = out.haplotaggedBamN
	haplotaggedBamidxN     = out.haplotaggedBamidxN
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
    haplotaggedBam         { path 'haplotagged_bam/tumor' }
    haplotaggedBamidx      { path 'haplotagged_bam/tumor' }
    haplotaggedBamN        { path 'haplotagged_bam/normal' }
    haplotaggedBamidxN     { path 'haplotagged_bam/normal' }
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


