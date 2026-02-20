# Lumos: Long-Read Somatic Variant Calling

This Nextflow (v25.10.0) workflow implements a full long-read somatic variant analysis pipeline.  
It supports tumor–normal and tumor-only configurations, running on [Biowulf](https://hpc.nih.gov/) (Slurm scheduler).  

Currently, the pipeline includes the following steps:

- **Alignment** with minimap2  
- **Small variant calling** with Clair3  
- **Phasing** with longphase  
- **Somatic SV calling** with Severus  
- **Copy-number analysis (CNA)** with Wakhan  
- **Somatic SNV calling** with DeepSomatic  
- **Methylation analysis** with modkit  

---

## 🚀 Quickstart

### Tumor–Normal Run

```bash
module load nextflow/25.10.0
nextflow run tumorNormalONT.nf   --reads_tumor tumor.bam   --reads_normal normal.bam   --reference hg38.fasta   --vntr vntr.bed   --clair3_model clair3_models/ont --cpgs cpgs.bed
```

### Tumor-Only Run

```bash
module load nextflow/25.10.0
nextflow run tumorOnlyONT.nf   --reads_tumor tumor.bam   --reference hg38.fasta   --vntr vntr.bed   --sv_pon PoN_1000G_hg38.tsv.gz   --clair3_model clair3_models/ont --cpgs cpgs.bed
```

> **Tip:** Always run inside an interactive session or with `sbatch` — **not** on the Biowulf head node.  

---

## 📥 Required Inputs

### Tumor–Normal Only

```
--reads_normal  Path to normal reads (unaligned bam)
--reads_tumor   Path to tumor reads (unaligned bam)
```

### Tumor-Only Only

```
--reads_tumor   Path to tumor reads (unaligned bam)
--sv_pon        Panel of Normals file (e.g. ./annot/PoN_1000G_hg38_extended.tsv.gz)
```

### For both modes

```
--reference     Reference FASTA  
--vntr          BED file of tandem repeats (must be ordered)(e.g. ./annot/human_GRCh38_no_alt_analysis_set.trf.bed)
--clair3_model  Path to Clair3 model  
--cpgs          CpG island BED file required in sv_cna_dmr and all modes (e.g. ./annot/hg38_cpg_cleaned.bed)
```

### Pre-existing bam alignment, use:
```
--aligned_input true --aligned_tumor BAM --aligned_tumor_bai BAI [--aligne_normal BAM --aligned_normal_bai BAI]
```

instead of `--reads_normal` / `--reads_tumor`

---

## ⚙️ Optional Parameters

```
--mode       sv_cna       Run only SV and CNA calling  
             sv_cna_dmr   Run SV, CNA, and DMR calling  
             all          Run SV, CNA, DMR, and somatic SNV calling (default)
--cosmic        Path to COSMIC genes in tsv format for Wakhan visualization
```


---

## ⚠️ Notes & Caveats

- **Input reads:** Must be in unmapped BAM format.  
  - Multiple files can be passed as a space-separated string wrapped in single quotes:  
    `--reads_tumor 'BAM1 BAM2'`  
  - Wildcards are supported:  
    `--reads_tumor 'BAM_DIR/*bam'`
- **Multiple runs:** Each Nextflow run must be launched from a unique working directory.  
  Nextflow creates a `work/` folder and `.nextflow.log` inside the run directory.
- **Resuming after failure:** Use the `-resume` flag (note single dash).  
  The workflow will attempt to reuse existing results.  
  Resume must be launched from the *same* working directory.  
- **Output directory change** By default, final outputs will be available in `lumos_out` directory.
  You can change it by adding `-output-dir XXX` argument (note single dash)
- **Nexflow version** Note that the workflow has been tested with Nextflow v25.10.0,
  and may not run with other versions due to recent synthax changes.

---

## 🛠️ Debugging & Development Tips

- To run locally, comment out  
  ```groovy
  process.executor = 'slurm'
  ```  
  in `nextflow.config`.
- Job status looks like:  
  ```
  [11/4852c3] tumorOnlyOntWorkflow:alignMinimap2 (1) [ 0%] 0 of 1
  ```
  Outputs are in the corresponding working directory under `work/11/4852c3...`.
- Useful files inside a job’s work directory:
  - `.command.sh` → exact command executed  
  - `.command.out` → stdout  
  - `.command.err` → stderr  
- Each process runs inside a Singularity/Docker container.  
  Custom container build scripts are available in the `docker/` folder.  

---

## 📚 Learning Resources

- New to Nextflow? Start with the excellent tutorial:  
  👉 https://training.nextflow.io/
  

**Notes**
- **Tumor-only:** omit `--reads_normal`; provide `--sv_pon` for Severus.
- **Modes:** 
  - `--mode sv_cna` runs **Severus** + **Wakhan**.  
  - `--mode sv_cna_dmr` adds **DMR** from modkit.  
  - `--mode all` runs all modules including **DeepSomatic**.
- **Consistency:** tumor and normal BAMs are haplotagged by the same method to ensure consistent phasing.
