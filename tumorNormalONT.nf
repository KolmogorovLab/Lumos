#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

include { callClair3; phaseLongphase; severusTumorNormal; wakhanCNATN; wakhanHapcorrectTN; modkitPileupAllele; deepsomaticTumorNormal } from "./processes/processes.nf"
include { alignMinimap2 as alignTumor } from "./processes/processes.nf"
include { alignMinimap2 as alignNormal } from "./processes/processes.nf"
include { haplotagWhatshap as haplotagNormal } from "./processes/processes.nf"
include { haplotagWhatshap as haplotagTumor } from "./processes/processes.nf"
include {modkitStats as modkitStats} from "./processes/processes.nf"
include {modkitStats as modkitStats2} from "./processes/processes.nf"
include {modkitStats as modkitStats3} from "./processes/processes.nf"
include {modkitPileup as modkitPileup} from "./processes/processes.nf"
include {modkitPileup as modkitPileupNormal} from "./processes/processes.nf"
include {modkitDMR as modkitDMRAllele} from "./processes/processes.nf"
include {modkitDMR as modkitDMR} from "./processes/processes.nf"

/*
 * Main workflow
 */
workflow tumorNormalOntWorkflow {
    take:
        readsTumor
        readsNormal
        reference
        vntrAnnotation
        clair3Model
	cpgs

    main:
        alignTumor(reference, readsTumor.collect())
        alignNormal(reference, readsNormal.collect())
        callClair3(alignNormal.out.bam, alignNormal.out.bam_idx, reference, alignNormal.out.ref_idx, clair3Model)
        phaseLongphase(alignNormal.out.bam, alignNormal.out.bam_idx, reference, 
                       alignNormal.out.ref_idx, callClair3.out.vcf)
	wakhanHapcorrectTN(alignMinimap2.out.bam, alignMinimap2.out.bam_idx, reference, phaseLongphase.out.phasedVcf)
        haplotagNormal(reference, alignNormal.out.ref_idx, wakhanHapcorrectTN.out.rephasedVcf, alignNormal.out.bam,
                       alignNormal.out.bam_idx)
        haplotagTumor(reference, alignTumor.out.ref_idx, wakhanHapcorrectTN.out.rephasedVcf, alignTumor.out.bam,
                      alignTumor.out.bam_idx)			   
        severusTumorNormal(haplotagTumor.out.bam, haplotagTumor.out.bam_idx,haplotagNormal.out.bam, haplotagNormal.out.bam_idx,
                           wakhanHapcorrectTN.out.rephasedVcf, vntrAnnotation)
        wakhanCNATN(alignMinimap2.out.bam, alignMinimap2.out.bam_idx, reference, wakhanHapcorrectTN.out.rephasedVcf, severusTumorNormal.out.severusSomaticVcf, wakhanHapcorrectTN.out.wakhanOutput)
        modkitPileupAllele(haplotagTumor.out.bam, haplotagTumor.out.bam_idx, reference, alignMinimap2.out.ref_idx)
        modkitPileupTumor(haplotagTumor.out.bam, haplotagTumor.out.bam_idx,reference, alignMinimap2.out.ref_idx)
	modkitPileupNormal(haplotagNormal.out.bam, haplotagNormal.out.bam_idx,reference, alignMinimap2.out.ref_idx)
        modkitDMRAllele(modkitPileupAllele.out.HP1bed, modkitPileupAllele.out.HP1bed_idx, modkitPileupAllele.out.HP2bed, modkitPileupAllele.out.HP2bed_idx, reference, alignNormal.out.ref_idx, cpgs)
	modkitDMR(modkitPileupNormal.out.pileupbed, modkitPileupNormal.out.pileupbed_idx, modkitPileup.out.pileupbed, modkitPileup.out.pileupbed_idx, reference, alignNormal.out.ref_idx, cpgs)
        modkitStats1(modkitPileupAllele.out.HP1bed, modkitPileupAllele.out.HP1bed_idx, reference, alignNormal.out.ref_idx, cpgs)
        modkitStats2(modkitPileupAllele.out.HP2bed, modkitPileupAllele.out.HP2bed_idx, reference, alignNormal.out.ref_idx, cpgs)
        modkitStats3(modkitPileup.out.pileupbed, modkitPileup.out.pileupbed_idx, reference, alignNormal.out.ref_idx, cpgs), reference, alignMinimap2.out.ref_idx)
        deepsomaticTumorNormal(alignTumor.out.bam, alignTumor.out.bam_idx, alignNormal.out.bam, alignNormal.out.bam_idx,
                               reference, alignTumor.out.ref_idx)

    emit:
        phasedVcf = phaseLongphase.out.phasedVcf
        haplotaggedTumor = haplotagTumor.out.bam
        haplotaggedTumorIdx = haplotagTumor.out.bam_idx
        haplotaggedNormal = haplotagNormal.out.bam
        haplotaggedNormalIdx = haplotagNormal.out.bam_idx
        severusFullOutput = severusTumorNormal.out.severusFullOutput
        wakhanFullOutput = wakhanTumorNormal.out.wakhanOutput
        modkitPileupAlleleBED1 = modkitPileupAllele.out.HP1bed
        modkitPileupAlleleBED2 = modkitPileupAllele.out.HP2bed
        modkitPileupOut = modkitPileup.out.pileupbed
	modkitPileupNormalOut = modkitPileupNormal.out.pileupbed
        modkitDMRAlleleOut = modkitDMRAllele.out.DMRbed
	modkitDMROut = modkitDMR.out.DMRbed
        modkitStatsOut = modkitStats.out.stats
        modkitStats2Out = modkitStats2.out.stats
        modkitStats3Out = modkitStats3.out.stats
	deepsomaticOutput = deepsomaticTumorNormal.out.deepsomaticOutput
	
    publish:
        phasedVcf >> "phased_vcf"
        haplotaggedTumor >> "haplotagged_bam_tumor"
        haplotaggedTumorIdx >> "haplotagged_bam_tumor"
        haplotaggedNormal >> "haplotagged_bam_normal"
        haplotaggedNormalIdx >> "haplotagged_bam_normal"
        severusFullOutput >> "severus"
        wakhanFullOutput >> "wakhan"
        deepsomaticOutput >> "deepsomatic"
	modkitPileupAlleleBED1 >> "methlation_allele"
        modkitPileupAlleleBED2 >> "methlation_allele"
	modkitDMRAlleleOut >> "methlation_allele"
        modkitPileupOut >> "methlation"
	modkitPileupNormalOut >> "methlation"
        modkitDMROut >> "methlation"
        modkitStatsOut >> "methlation"
        modkitStats2Out >> "methlation"
        modkitStats3Out >> "methlation"
}

/*
 * Entry point
 */
workflow {
    if (!params.reads_tumor || !params.reads_normal || !params.reference || !params.outdir || 
        !params.vntr || !params.clair3_model || !params.cpgs) {
        error """
              ERROR: Some required arguments are not defined.
              Usage: tumorNormalONT.nf --reads_tumor PATH --reads_normal PATH --reference PATH --outdir PATH 
                                       --vntr PATH --clair3_model PATH --alignment [true, false] --mode[all, sv_cna, sv_cna_dmr]
              """.stripIndent()
    }

    params.mode = params.mode ?: 'all'
    params.alignment = params.alignment ?: 'true'

    tumorChannel = Channel.fromPath(params.reads_tumor.split(" ").toList(), checkIfExists: true)
    tumorChannel.view{it -> "Tumor reads: $it"}

    normalChannel = Channel.fromPath(params.reads_normal.split(" ").toList(), checkIfExists: true)
    normalChannel.view{it -> "Normal reads: $it"}

    tumorNormalOntWorkflow(tumorChannel, normalChannel,
                           Channel.fromPath(params.reference, checkIfExists: true), 
                           Channel.fromPath(params.vntr, checkIfExists: true),
			   Channel.fromPath(params.cpgs, checkIfExists: true),
                           Channel.fromPath(params.clair3_model, checkIfExists: true))
}

output {
    directory params.outdir
    mode "copy"
}	
