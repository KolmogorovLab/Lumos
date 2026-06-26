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
        path reads, stageAs: "input_reads/*"
    output:
        path 'aligned.bam', emit: bam
        path 'aligned.bam.bai', emit: bam_idx
        path "${ref}.fai", emit: ref_idx
          
    script:
        """  
        samtools cat input_reads/* | \
          samtools fastq -TMm,Ml,MM,ML - | \
          minimap2 -ax map-ont -k 17 -t ${task.cpus} -K 1G -y --eqx ${ref} - | \
          samtools sort -@4 -m 4G > aligned.bam
        samtools index -@8 aligned.bam
        samtools faidx ${ref}
        """
}   


process fastaIndex {
  tag "${reference?.name}"
  container 'docker://quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1'
  cpus 1
  memory '4G'
  executor 'local'

  input:
    path reference
  output:
    path "${reference}.fai", emit: ref_idx
  script:
  """
  samtools faidx ${reference}
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
        path alignedBam, stageAs: "input.bam"
        path indexedBai, stageAs: "input.bam.bai"

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

def MODKIT_DOCKER = 'docker://mkolmogo/modkit:0.4.1'

process modkitDMR{
        container MODKIT_DOCKER
        cpus 28
        memory '64 G'
        time '10.h'

        input:
            path tumorBed
		    path tumorBedIdx
		    path normalBed
		    path normalBedIdx
            path reference
            path referenceIdx
            path cpgbed

        output:
                path 'dmr.bed', emit:DMRbed
        script:
            """
            modkit dmr pair -a ${tumorBed} -b ${normalBed} -o ./dmr.bed -r ${cpgbed} --ref ${reference} --base C --threads ${task.cpus}
            """
}

process modkitStats{
        container MODKIT_DOCKER
        cpus 28
        memory '64 G'
        time '10.h'

        input:
            path tumorBed
		    path tumorBedIdx
            path reference
            path referenceIdx
            path cpgbed

        output:
        	path "${tumorBed.simpleName}.stats.tsv", emit: stats

        script:
            """
            modkit stats ${tumorBed} --regions ${cpgbed} -o ./${tumorBed.simpleName}.stats.tsv --mod-codes "h,m" --threads ${task.cpus}
            """
}

process modkitPileupAllele{
        container MODKIT_DOCKER
        cpus 28
        memory '64 G'
        time '10.h'

        input:
            path tumorBam
            path tumorBamIdx
            path reference
            path referenceIdx

	    output:
            path 'haplotagged_1.bed.gz', emit:HP1bed
            path 'haplotagged_2.bed.gz', emit:HP2bed
            path 'haplotagged_1.bed.gz.tbi', emit:HP1bed_idx
            path 'haplotagged_2.bed.gz.tbi', emit:HP2bed_idx

        script:
            """
            modkit pileup ${tumorBam} . -t ${task.cpus} --combine-strands --cpg --ref ${reference} --no-filtering --partition-tag HP --prefix haplotagged
            bgzip haplotagged_1.bed
            tabix -p bed haplotagged_1.bed.gz
            bgzip haplotagged_2.bed
            tabix -p bed haplotagged_2.bed.gz
            """     
}


process modkitPileup{
        container MODKIT_DOCKER
        cpus 28
        memory '64 G'
        time '10.h'

        input:
            path tumorBam
            path tumorBamIdx
            path reference
            path referenceIdx
		    path cpgbed

	    output:
            path 'pileup.bed.gz', emit:pileupbed
            path 'pileup.bed.gz.tbi', emit:pileupbed_idx
            path 'pileup_cpg_subset.bed.gz',emit:pileupbed_subset
            path 'pileup_cpg_subset.bed.gz.tbi', emit:pileupbed_subset_idx

        script:
            """
            modkit pileup ${tumorBam} pileup.bed -t ${task.cpus} --combine-strands --cpg --ref ${reference} --no-filtering
            bgzip pileup.bed
            tabix -p bed pileup.bed.gz
            bedtools intersect -a pileup.bed.gz -b ${cpgbed} -wa -wb > pileup_cpg.bed 
            cut -f1-5,12,19-21 pileup_cpg.bed > pileup_cpg_subset.bed
            bgzip pileup_cpg_subset.bed
            tabix -p bed pileup_cpg_subset.bed.gz
            """     
}

def SEVERUS_DOCKER = 'docker://mkolmogo/severus:1.8dev'

process severusTumorOnly {
    container SEVERUS_DOCKER
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
        path 'severus_out/*', arity: '3..*', emit: severusFullOutput
        path 'severus_out/somatic_SVs/severus_somatic.vcf', emit: severusSomaticVcf

    script:
        """
        tabix ${phasedVcf}
        severus --target-bam ${tumorBam} --out-dir severus_out -t ${task.cpus} --phasing-vcf ${phasedVcf} \
            --vntr-bed ${vntrBed} --PON ${panelOfNormals} --output-read-ids --min-reference-flank 0 --single-bp --resolve-overlaps --max-unmapped-seq 7000 --between-junction-ins 
        """
}

process severusTumorNormal {
    container SEVERUS_DOCKER
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

def WAKHAN_DOCKER = 'docker://mkolmogo/wakhan:0.4.3'
def WAKHAN_BIN = 50000

process wakhanHapcorrect {
    def genomeName = "Sample"

    container WAKHAN_DOCKER
    cpus 16
    memory '64 G'
    time '14.h'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path reference
        path tumorSmallPhasedVcf

    output:
        path 'wakhan_hapcorrect/*', arity: '3..*', emit: wakhanHPOutput
        path 'wakhan_hapcorrect/phasing_output/rephased.vcf.gz', emit: rephasedVcf

    script:
        """
        tabix ${tumorSmallPhasedVcf}
        wakhan hapcorrect --threads ${task.cpus} --reference ${reference} --target-bam ${tumorBam} --tumor-phased-vcf ${tumorSmallPhasedVcf} \
          --genome-name Sample --out-dir-plots wakhan_hapcorrect --bin-size ${WAKHAN_BIN}  --phaseblocks-enable \
          --contigs ${params.contigs ?: 'chr1-22,chrX'} --copynumbers-subclonal-enable
        """
}

process wakhanCNA {
    def genomeName = "Sample"
    container WAKHAN_DOCKER
    cpus 16
    memory '64 G'
    time '14.h'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path reference
        path tumorSmallPhasedVcf
        path severusSomaticVcf
	    path hapcorrect_out, stageAs: "hapcorrect_input/*"
        path cosmic

    output:
        path "wakhan_cna", emit: wakhanOutput

    script:
        def cosmic_arg = cosmic.name != 'NOFILE' ? "--cancer-genes ${cosmic}" : ""
        """
        tabix ${tumorSmallPhasedVcf}
        cp -rL hapcorrect_input wakhan_cna
        wakhan cna --threads ${task.cpus} --reference ${reference} --target-bam ${tumorBam} --tumor-phased-vcf ${tumorSmallPhasedVcf} \
          --genome-name Sample --use-sv-haplotypes --out-dir-plots wakhan_cna --bin-size ${WAKHAN_BIN}  --breakpoints ${severusSomaticVcf} --phaseblocks-enable \
          --contigs ${params.contigs ?: 'chr1-22,chrX'} --copynumbers-subclonal-enable ${cosmic_arg}
	"""
}

process wakhanHapcorrectTN {
    def genomeName = "Sample"

    container WAKHAN_DOCKER
    cpus 16
    memory '64 G'
    time '14.h'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path reference
        path tumorSmallPhasedVcf

    output:
        path 'wakhan_hapcorrect/*', arity: '3..*', emit: wakhanHPOutput
        path 'wakhan_hapcorrect/phasing_output/rephased.vcf.gz', emit: rephasedVcf

    script:
        """
        tabix ${tumorSmallPhasedVcf}
        wakhan hapcorrect --threads ${task.cpus} --reference ${reference} --target-bam ${tumorBam} --normal-phased-vcf ${tumorSmallPhasedVcf} \
          --genome-name Sample --out-dir-plots wakhan_hapcorrect --bin-size ${WAKHAN_BIN}  --phaseblocks-enable \
          --contigs ${params.contigs ?: 'chr1-22,chrX'} --copynumbers-subclonal-enable
        """
}

process wakhanCNATN {
    def genomeName = "Sample"
    container WAKHAN_DOCKER
    cpus 16
    memory '64 G'
    time '14.h'

    input:
        path tumorBam, stageAs: "tumor.bam"
        path tumorBamIdx, stageAs: "tumor.bam.bai"
        path reference
        path tumorSmallPhasedVcf
        path severusSomaticVcf
        path hapcorrect_out, stageAs: "hapcorrect_input/*"
        path cosmic

    output:
        path "wakhan_cna", emit: wakhanOutput

    script:
        def cosmic_arg = cosmic.name != 'NOFILE' ? "--cancer-genes ${cosmic}" : ""
        """
        tabix ${tumorSmallPhasedVcf}
        cp -rL hapcorrect_input wakhan_cna
        wakhan cna --threads ${task.cpus} --reference ${reference} --target-bam ${tumorBam} --normal-phased-vcf ${tumorSmallPhasedVcf} \
          --use-sv-haplotypes --genome-name Sample --out-dir-plots wakhan_cna --bin-size ${WAKHAN_BIN}  --breakpoints ${severusSomaticVcf} --phaseblocks-enable \
          --contigs ${params.contigs ?: 'chr1-22,chrX'} --copynumbers-subclonal-enable ${cosmic_arg}
	"""
}

def DEEPSOMATIC_DOCKER = 'docker://google/deepsomatic:1.9.0'

process deepsomaticTumorOnly {
    def genomeName = "Sample"
    def outDir = "deepsomatic_out"

    container DEEPSOMATIC_DOCKER
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
    container DEEPSOMATIC_DOCKER
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
