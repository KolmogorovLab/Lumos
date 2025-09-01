# Lumos: Long-Read Somatic Variant Calling

This Nextflow (v25.04.2) workflow implements a full long-read somatic variant analysis pipeline.  
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
module load nextflow
nextflow run tumorNormalONT.nf   --reads_tumor tumor.bam   --reads_normal normal.bam   --reference hg38.fasta   --outdir results   --vntr vntr.bed   --clair3_model clair3_models/ont
```

### Tumor-Only Run

```bash
module load nextflow
nextflow run tumorOnlyONT.nf   --reads tumor.bam   --reference hg38.fasta   --outdir results   --vntr vntr.bed   --sv_pon PoN_1000G_hg38.tsv.gz   --clair3_model clair3_models/ont
```

> **Tip:** Always run inside an interactive session or with `sbatch` — **not** on the Biowulf head node.  

---

## 📥 Required Inputs

### Common

```
--reads_tumor   Path to tumor BAM file(s), must be indexed  
--outdir        Output directory  
--reference     Reference FASTA  
--vntr          BED file of tandem repeats (must be ordered)(e.g. ./annot/human_GRCh38_no_alt_analysis_set.trf.bed)
--clair3_model  Path to Clair3 model  
--cpgs          CpG island BED file (e.g. ./annot/hg38_cpg_cleaned.bed)
```

### Tumor–Normal Only

```
--reads_normal  Path to normal BAM file (must be indexed)
```

### Tumor-Only Only

```
--sv_pon        Panel of Normals file (e.g. ./annot/PoN_1000G_hg38_extended.tsv.gz)
```

---

## ⚙️ Optional Parameters

```
--mode       sv_cna       Run only SV and CNA calling  
             sv_cna_dmr   Run SV, CNA, and DMR calling  
             all          Run SV, CNA, DMR, and somatic SNV calling (default)

--aligned    Provide pre-aligned BAM files. Use:
             --tumor_bam, --tumor_bai
             --normal_bam, --normal_bai
             instead of --reads_tumor / --reads_normal
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
- **Resuming after failure:** Use the `-resume` flag (single dash).  
  The workflow will attempt to reuse existing results.  
  Resume must be launched from the *same* working directory.  

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
  
## 🗺️ Workflow Diagram (Mermaid)

```mermaid
flowchart TD
  %% ============ Inputs ============
  subgraph I[Inputs]
    TUMOR[Tumor reads/BAM (indexed)]
    NORMAL[Normal reads/BAM (indexed, optional)]
    REF[Reference FASTA]
    VNTR[VNTR BED]
    CPG[CpG islands BED]
    PON[SV Panel of Normals (tumor-only)]
    CLAIR3[Clair3 model]
  end

  %% ============ Core Steps ============
  ALIGN[Alignment: minimap2]
  T_BAM[Aligned Tumor BAM]
  N_BAM[Aligned Normal BAM]
  SMALLVAR[Small variants: Clair3]
  PHASE[Phasing: longphase]
  HAPTAG[Haplotagging (tag tumor/normal consistently)]
  SV[Somatic SV calling: Severus]
  CNA[CNA calling: Wakhan]
  METH[Methylation: modkit]
  DMR[DMR calling (modkit dmr/stats)]
  SNV[Somatic SNV calling: DeepSomatic]
  REPORTS[QC & Reports]
  OUT[Final Outputs]

  %% ============ Flow ============
  TUMOR --> ALIGN
  NORMAL --> ALIGN
  REF --> ALIGN
  ALIGN --> T_BAM
  ALIGN --> N_BAM

  T_BAM --> SMALLVAR
  N_BAM --> SMALLVAR
  CLAIR3 --> SMALLVAR
  SMALLVAR --> PHASE
  PHASE --> HAPTAG

  HAPTAG --> SV
  HAPTAG --> CNA
  HAPTAG --> SNV
  HAPTAG --> METH

  VNTR -. annotation .- SV
  VNTR -. annotation .- SMALLVAR
  CPG -. regions .- METH

  METH --> DMR

  %% Tumor-only branch hint
  PON -. used in tumor-only .-> SV

  %% Outputs
  SV --> OUT
  CNA --> OUT
  SNV --> OUT
  DMR --> OUT
  OUT --> REPORTS
```

**Notes**
- **Tumor-only:** omit `--reads_normal`; provide `--sv_pon` for Severus.
- **Modes:** 
  - `--mode sv_cna` runs **Severus** + **Wakhan**.  
  - `--mode sv_cna_dmr` adds **DMR** from modkit.  
  - `--mode all` runs all modules including **DeepSomatic**.
- **Consistency:** tumor and normal BAMs are haplotagged by the same method to ensure consistent phasing.