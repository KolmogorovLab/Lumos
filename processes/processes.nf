#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Align reads using minimap2, sort BAM using samtools, and create BAM index
 */
process alignMinimap2 {

    container 'docker://quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1'
    cpus 28
    memory '128 GB'
    time '24.h'

    input:
        path ref
        path reads

    output:
        path 'aligned.bam', emit: bam
        path 'aligned.bam.bai', emit: bam_idx
        path "${ref}.fai", emit: ref_idx
          
    script:
        """  
        samtools cat ${reads} | \
          samtools fastq -TMm,Ml,MM,ML - | \
          minimap2 -ax map-ont -k 17 -t ${task.cpus} -K 1G -y --eqx ${ref} - | \
          samtools sort -@4 -m 4G > aligned.bam
        samtools index -@8 aligned.bam
        samtools faidx ${ref}
        """
}   
    
process callClair3 {
    container 'docker://hkubal/clair3:v1.0.11'
    cpus 28
    memory '128 G'
    time '24.h'

    input:
        path alignedBam
        path indexedBai
        path reference
        path referenceIdx
        path modelPath

    output:
        path 'clair3_output/merge_output.vcf.gz', emit: vcf

    
    script:
        """
        /opt/bin/run_clair3.sh \
            --bam_fn=${alignedBam} \
            --ref_fn=${reference} \
            --threads=${task.cpus} \
            --platform="ont" \
            --model_path=${modelPath} \
            --output="clair3_output" \
        """
}

process phaseLongphase {

    container 'docker://mkolmogo/longphase:1.7.3'
    cpus 10
    memory '64 G'
    time '4.h'

    input:
        path alignedBam
        path indexedBai
        path reference
        path referenceIdx
        path vcf

    output:
        path 'longphase.vcf.gz', emit: phasedVcf

    
    script:
        """
        longphase phase -s ${vcf} -b ${alignedBam} -r ${reference} -t ${task.cpus} -o longphase --ont
        bgzip longphase.vcf
        """
}

process haplotagWhatshap {
    container 'docker://mkolmogo/whatshap:2.3'
    cpus 8
    memory '64 G'
    time '10.h'

    input:
        path reference
        path referenceIdx
        path phasedVcf
        path alignedBam
        path indexedBai

    output:
        path 'haplotagged.bam', emit: bam
        path 'haplotagged.bam.bai', emit: bam_idx

	    script:
        """
        tabix ${phasedVcf}
        whatshap haplotag --reference ${reference} ${phasedVcf} ${alignedBam} -o 'haplotagged.bam' --ignore-read-groups \
            --tag-supplementary --skip-missing-contigs --output-threads 4
        samtools index -@8 haplotagged.bam
        """
}

process severusTumorOnly {
    container 'docker://gokcekeskus/severus:dev_134407c'
    cpus 28
    memory '128 G'
    time '8.h'

    input:
        path tumorBam
        path tumorBamIdx
        path phasedVcf
        path vntrBed
        path panelOfNormals

    output:
        path 'severus_out_15/*', arity: '3..*', emit: severusFullOutput
        path 'severus_out_15/somatic_SVs/severus_somatic.vcf', emit: severusSomaticVcf

    script:
        """
        tabix ${phasedVcf}
        severus --target-bam ${tumorBam} --out-dir severus_out_15 -t ${task.cpus} --phasing-vcf ${phasedVcf} \
            --vntr-bed ${vntrBed} --PON ${panelOfNormals} --output-read-ids --min-reference-flank 0 --single-bp --resolve-overlaps --max-unmapped-seq 7000 --between-junction-ins 
        """
}

process severusTumorNormal {
    container 'docker://gokcekeskus/severus:dev_134407c'
    cpus 28
    memory '128 G'
    time '8.h'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path normalBam, stageAs: "normal.bam"
        path normalBamIdx, stageAs: "normal.bam.bai"
        path phasedVcf
        path vntrBed

    output:
        path 'severus_out/*', arity: '3..*', emit: severusFullOutput
        path 'severus_out/somatic_SVs/severus_somatic.vcf', emit: severusSomaticVcf

    script:
        """
        tabix ${phasedVcf}
        severus --target-bam ${tumorBam} --control-bam ${normalBam} --out-dir severus_out -t ${task.cpus} --phasing-vcf ${phasedVcf} \
            --vntr-bed ${vntrBed} --single-bp --resolve-overlaps --max-unmapped-seq 7000 --between-junction-ins 
        """
}

process wakhanTumorOnly {
    def genomeName = "Sample"

    container 'docker://gokcekeskus/wakhan:dev_3a3ea29'
    cpus 16
    memory '64 G'
    time '14.h'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path reference
        path tumorSmallPhasedVcf
        path severusSomaticVcf
	path gene_list

    output:
        path 'wakhan_out_88e88ca_npcpd/*', arity: '3..*', emit: wakhanOutput
        path 'wakhan_out_88e88ca_nocpd/phasing_output/Sample.rephased.vcf.gz', emit: rephasedVcf

    script:
        """
        tabix ${tumorSmallPhasedVcf}
        wakhan --threads ${task.cpus} --reference ${reference} --target-bam ${tumorBam} --tumor-vcf ${tumorSmallPhasedVcf} \
          --genome-name Sample --out-dir-plots wakhan_out_88e88ca_nocpd --bin-size 10000 --breakpoints ${severusSomaticVcf} \
		 --ploidy-range 1-6 --loh-enable --contigs chr1-22,chrX --user-input-genes ${gene_list} --copynumbers-subclonal-enable 
        """
}

process wakhanTumorNormal {
    def genomeName = "Sample"

    container 'docker://gokcekeskus/wakhan:dev_3a3ea29'
    cpus 16
    memory '64 G'
    time '14.h'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path reference
        path normalSmallPhasedVcf
        path severusSomaticVcf

    output:
        path 'wakhan_out/*', arity: '3..*', emit: wakhanOutput

    script:
        """
        tabix ${normalSmallPhasedVcf}
        wakhan --threads ${task.cpus} --reference ${reference} --target-bam ${tumorBam} --normal-phased-vcf ${normalSmallPhasedVcf} \
          --genome-name Sample --out-dir-plots wakhan_out --breakpoints ${severusSomaticVcf} --bin-size-snps 100000 --bin-size 10000 \
		--ploidy-range 1-6 --cpd-internal-segments --copynumbers-subclonal-enable
        """
}

process deepsomaticTumorOnly {
    def genomeName = "Sample"
    def outDir = "deepsomatic_out"

    container 'docker://google/deepsomatic:1.7.0'
    cpus 56
    memory '240 G'
    time '48.h'
    clusterOptions '--exclusive'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path reference
        path referenceIdx

    output:
        path 'deepsomatic_out/ds.merged.vcf.gz', emit: deepsomaticOutput

    script:
        """
        ds_parallel_tumor_only.sh ${tumorBam} ${reference} ${outDir} ${genomeName}
        """
}

process deepsomaticTumorNormal {
    container 'docker://google/deepsomatic:1.7.0'
    cpus 56
    memory '240 G'
    time '48.h'
    clusterOptions '--exclusive'

    def genomeName = "Sample"
    def outDir = "deepsomatic_out"

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path normalBam, stageAs: "normal.bam"
        path normalBamIdx, stageAs: "normal.bam.bai"
        path reference
        path referenceIdx

    output:
        path 'deepsomatic_out/ds.merged.vcf.gz', emit: deepsomaticOutput

    script:
        """
        ds_parallel_tumor_normal.sh ${tumorBam} ${normalBam} ${reference} ${outDir} ${genomeName}-T ${genomeName}-N
        """
}
