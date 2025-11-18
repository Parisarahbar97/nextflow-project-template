Demux + QC Intersection Pipeline (dsl2)
=======================================

This project runs an end-to-end demultiplexing workflow for all pools, producing high-confidence singlet labels as the intersection of:
- popscle demuxlet (harmonized 1000G pileup + donor VCF, alpha 0/0.5, prior 0.05, ERR 0.15)
- scQC-flow transcriptomic QC (filters: `doublet_class == "singlet"`, `nuclear_fraction >= 0.4`, `percent.mt <= 20`)

Inputs (pools.csv)
------------------
`pool,bam,barcodes,donor_vcf,qc_rds,outdir` for each pool (S1A–S19B). Example paths used here:
- BAM: `/home/pr422/RDS/live/Users/Parisa/alex_output/epilep_cellranger_outputs/${POOL}_mapped/outs/possorted_genome_bam.bam`
- Barcodes (CellBender whitelist): `/home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/qc/output_latest/${POOL}/${POOL}_cellbender_output/cellbender_out_cell_barcodes.csv`
- Donor VCF: `/home/pr422/RDS/live/Users/Parisa/alex_output/genotype_inputs2/${POOL}.vcf.gz`
- QC Seurat RDS: `/home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/qc/output_latest/${POOL}/${POOL}_seurat_object.rds`
- Outdir: `/home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/demuxlet_results/${POOL}`

Global defaults (main.nf)
-------------------------
- Pop VCF: `/home/pr422/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.gz`
- FASTA: `/home/pr422/RDS/live/Users/Parisa/refrence_genome/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa`
- Demuxlet: `alpha 0,0.5`, `--doublet-prior 0.05`, `--geno-error-offset 0.15`, `--geno-error-coeff 0.0`

Pipeline steps (per pool)
-------------------------
1) Normalize pop VCF (MAF≥0.05, biallelic SNVs)  
2) Normalize donor VCF (bcftools +fixref/norm, biallelic SNVs)  
3) Intersect pop/donor → `sites.intersect.vcf.gz`  
4) Reheader to BAM order → `sites.intersect.bamorder.vcf.gz`  
5) Harmonize donor to sites → `donor.harmonized.to_sites.vcf.gz`  
6) Pileup: `scripts/01_run_pileup_1000G.sh` (CB/UB, MQ20, BQ13, cap BQ40)  
7) Demuxlet: `scripts/02_run_demuxlet_sweep.sh` (alpha 0/0.5, prior 0.05, ERR 0.15) using the barcode whitelist  
8) QC ∩ demux intersection: `scripts/qc_intersect.R` → `<pool>_final_labels.tsv` copied to outdir  

Containers
----------
- QC_INTERSECT runs in `ghcr.io/johnsonlab-ic/sc_analysis:latest` (set in `nextflow.config`, docker enabled for local profile).  
- Pileup/demuxlet call their own Docker images inside the bash scripts (`parisa/demux:2.1`).

How to run
----------
```bash
cd /home/pr422/RDS/live/Users/Parisa/demux-qc-flow
nextflow run main.nf -profile local   --pools_csv pools.csv   # local profile (Docker on)
# or
nextflow run main.nf -profile standard --pools_csv pools.csv   # alias of local
```
Outputs per pool land in the `outdir` from `pools.csv` (e.g., `/home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/demuxlet_results/<POOL>`), including `<pool>_final_labels.tsv`.
