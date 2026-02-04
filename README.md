# ChIP-seq Analysis Pipeline

A reproducible Nextflow pipeline for ChIP-seq data analysis, encompassing quality control, read alignment, peak calling, and signal visualization using containerized tools.

## Overview

This pipeline performs comprehensive ChIP-seq analysis from raw FASTQ files to reproducible peak identification and signal quantification. Built with Nextflow and containerized tools via Singularity, it ensures complete reproducibility across different computing environments.

## Pipeline Workflow

The pipeline consists of the following major steps:

1. **Quality Control** - Assessment of raw sequencing data
2. **Read Trimming** - Adapter removal and quality filtering
3. **Read Alignment** - Alignment to reference genome
4. **Post-alignment Processing** - Sorting, indexing, and filtering
5. **Peak Calling** - Identification of enriched regions
6. **Reproducible Peak Identification** - Cross-replicate peak validation
7. **Peak Annotation** - Genomic feature assignment and motif discovery
8. **Signal Visualization** - Coverage track generation and correlation analysis

## Requirements

### Software Dependencies

All tools are provided through containerized environments (`ghcr.io/bf528`) executed via Singularity:

- [Nextflow](https://www.nextflow.io/) (≥21.0)
- [Singularity](https://sylabs.io/singularity/) (≥3.7)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.12.1)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (v0.39)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v2.5.4)
- [SAMtools](http://www.htslib.org/) (v1.21)
- [MultiQC](https://multiqc.info/) (v1.25)
- [HOMER](http://homer.ucsd.edu/homer/) (v5.1)
- [BEDtools](https://bedtools.readthedocs.io/) (v2.31.1)
- [deepTools](https://deeptools.readthedocs.io/) (v3.5.5)

### Reference Data Requirements

- Human reference genome: GRCh38 (gencode v47)
- Gene annotations: UCSC Genome Browser RefSeq (refGene)
- Blacklist regions: hg38 blacklist v2 (Boyle Lab)

### Input Requirements

- Single-end or paired-end FASTQ files
- At least two biological replicates per condition (recommended)
- Corresponding input/control samples

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/chipseq-pipeline.git
cd chipseq-pipeline

# Install Nextflow (if not already installed)
curl -s https://get.nextflow.io | bash

# Install Singularity (if not already installed)
# See: https://sylabs.io/guides/3.0/user-guide/installation.html
```

## Usage

### Basic Usage

```bash
nextflow run main.nf \
  --input samples.csv \
  --genome GRCh38 \
  --outdir results \
  -profile singularity
```

### Input Format

Create a CSV file (`samples.csv`) with the following format:

```csv
sample,fastq,replicate,type
sample1_rep1,/path/to/sample1_rep1.fastq.gz,1,treatment
sample1_rep2,/path/to/sample1_rep2.fastq.gz,2,treatment
input_rep1,/path/to/input_rep1.fastq.gz,1,control
input_rep2,/path/to/input_rep2.fastq.gz,2,control
```

For paired-end data:
```csv
sample,fastq_1,fastq_2,replicate,type
sample1_rep1,/path/to/sample1_rep1_R1.fastq.gz,/path/to/sample1_rep1_R2.fastq.gz,1,treatment
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to input CSV file | Required |
| `--genome` | Reference genome (GRCh38, GRCh37) | `GRCh38` |
| `--outdir` | Output directory | `./results` |
| `--blacklist` | Path to blacklist BED file | Auto-detected |
| `--peak_overlap` | Minimum reciprocal overlap for reproducible peaks | `0.5` |
| `--flank_size` | Flanking region size for TSS/TTS plots (bp) | `2000` |

### Advanced Options

```bash
nextflow run main.nf \
  --input samples.csv \
  --genome GRCh38 \
  --peak_overlap 0.6 \
  --flank_size 3000 \
  --outdir results \
  -profile singularity \
  -resume
```

## Pipeline Details

### 1. Quality Control (FastQC)

Raw FASTQ files are assessed for:
- Per-base sequence quality
- Adapter content
- GC content distribution
- Sequence duplication levels
- Overrepresented sequences

All analyses use default parameters.

### 2. Read Trimming (Trimmomatic)

Adapter sequences and low-quality bases are removed using default parameters:
- Adapter removal
- Quality trimming
- Minimum length filtering

### 3. Alignment (Bowtie2)

Reads are aligned to GRCh38 (gencode v47) using Bowtie2 with default parameters optimized for ChIP-seq data.

```bash
bowtie2 -x genome_index -U input.fastq.gz -S output.sam
```

### 4. Post-alignment Processing (SAMtools)

Aligned reads are processed for downstream analysis:
- **Sorting**: Coordinate-based sorting
- **Indexing**: BAM file indexing
- **Filtering**: Removal of unmapped and low-quality reads
- **Statistics**: flagstat utility for alignment metrics

All operations use default parameters.

### 5. Quality Report Aggregation (MultiQC)

FastQC, Trimmomatic, Bowtie2, and SAMtools metrics are aggregated into a single comprehensive HTML report.

### 6. Peak Calling (HOMER)

Enriched regions are identified using HOMER with default parameters:

```bash
findPeaks tags/ -i control_tags/ -style factor -o peaks.txt
```

### 7. Reproducible Peak Identification

Peaks are validated across biological replicates:

1. **Replicate Intersection**: Peaks from both replicates are intersected using BEDtools
2. **Overlap Threshold**: ≥50% reciprocal overlap required
   - Stringent enough to avoid false positives
   - Permissive enough to capture true biological reproducibility
3. **Blacklist Filtering**: Peaks overlapping hg38 blacklist v2 regions are removed

```bash
bedtools intersect -a rep1_peaks.bed -b rep2_peaks.bed -f 0.5 -r > reproducible_peaks.bed
bedtools intersect -a reproducible_peaks.bed -b blacklist.bed -v > final_peaks.bed
```

### 8. Peak Annotation and Motif Discovery (HOMER)

- Genomic feature annotation (promoter, intron, exon, intergenic)
- Distance to nearest TSS
- *De novo* motif discovery
- Known motif enrichment analysis

### 9. Signal Visualization (deepTools)

#### Coverage Tracks
BigWig files are generated from sorted BAM files using bamCoverage with default parameters:

```bash
bamCoverage -b input.bam -o output.bw
```

#### Signal Profiles
Coverage is computed around TSS and TTS using computeMatrix:
- **Mode**: scale-regions
- **Flanking regions**: ±2000 bp (customizable)
- **Reference**: UCSC RefSeq (refGene) annotations

#### Sample Correlation
Pearson correlation analysis is performed:
- **Tool**: multiBigWigSummary and plotCorrelation
- **Output**: Correlation heatmaps and scatterplots

All deepTools analyses use default parameters.

## Output Structure

```
results/
├── fastqc/                      # Raw read quality reports
│   ├── sample1_fastqc.html
│   └── sample1_fastqc.zip
├── trimmed/                     # Trimmed FASTQ files
│   └── sample1_trimmed.fastq.gz
├── aligned/                     # Alignment files
│   ├── sample1.bam
│   ├── sample1.bam.bai
│   └── sample1.flagstat
├── multiqc/                     # Aggregated QC report
│   └── multiqc_report.html
├── peaks/                       # Peak calling results
│   ├── sample1_peaks.txt
│   ├── sample1_peaks.bed
│   └── reproducible_peaks.bed
├── annotations/                 # Peak annotations
│   ├── peak_annotations.txt
│   ├── motif_results/
│   └── known_motifs/
├── bigwig/                      # Coverage tracks
│   └── sample1.bw
└── visualization/               # Plots and matrices
    ├── tss_signal_matrix.gz
    ├── tts_signal_matrix.gz
    ├── signal_profile.pdf
    ├── correlation_heatmap.pdf
    └── correlation_scatter.pdf
```

## Example Dataset

This pipeline was developed and validated using:
- **Data source**: NCBI GEO (GSE75070)
- **Organism**: Human (*Homo sapiens*)
- **Genome**: GRCh38 (gencode v47)
- **Design**: Biological replicates with matched input controls

## Reproducibility

All analyses are performed using containerized tools from `ghcr.io/bf528` executed through Singularity, ensuring:
- Consistent software versions across environments
- Portable execution on HPC clusters and cloud platforms
- Complete reproducibility of results

## Best Practices

### Quality Control Thresholds
- Minimum read quality: Q20
- Adapter contamination: <5%
- Library complexity: Check duplication rates

### Peak Calling Recommendations
- Always use matched input/control samples
- Minimum 2 biological replicates recommended
- Consider using IDR (Irreproducible Discovery Rate) for additional validation

### Replicate Overlap Threshold
The default 50% reciprocal overlap threshold balances:
- **Stringency**: Reduces false positive peaks
- **Sensitivity**: Captures genuine biological reproducibility
- Adjust based on experimental design and replicate concordance


### Tool Citations

- **FastQC**: Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data.
- **Trimmomatic**: Bolger et al. (2014) Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics.
- **Bowtie2**: Langmead & Salzberg (2012) Fast gapped-read alignment with Bowtie 2. Nature Methods.
- **SAMtools**: Danecek et al. (2021) Twelve years of SAMtools and BCFtools. GigaScience.
- **MultiQC**: Ewels et al. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics.
- **HOMER**: Heinz et al. (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Molecular Cell.
- **BEDtools**: Quinlan & Hall (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics.
- **deepTools**: Ramírez et al. (2016) deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research.

### Data Source

Data analyzed with this pipeline:
- **GEO Accession**: GSE75070
- **NCBI GEO**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75070

## Troubleshooting

### Common Issues

**Issue**: Singularity container fails to pull
```bash
# Solution: Manually pull the container
singularity pull docker://ghcr.io/bf528/chipseq-tools:latest
```

**Issue**: Memory errors during peak calling
```bash
# Solution: Increase memory in nextflow.config
process {
    withName: 'HOMER_FINDPEAKS' {
        memory = '32 GB'
    }
}
```

**Issue**: Low replicate concordance
- Check QC metrics for failed samples
- Verify input controls are appropriate
- Consider biological variability
- Adjust peak_overlap threshold if needed

## References

- ENCODE Project Consortium. (2012) An integrated encyclopedia of DNA elements in the human genome. Nature.
- Landt et al. (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Research.
- Amemiya et al. (2019) The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Scientific Reports.
