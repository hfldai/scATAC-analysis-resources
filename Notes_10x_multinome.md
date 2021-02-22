# 10x multinome of ATAC-seq and RNA-seq - data format & downstream analysis

## 10x Single-library analysis pipeline (software: `cellranger-arc`)
<img src="https://support.10xgenomics.com/img/cellranger-arc/algo_chart.png" width=600x>

## Raw data formats
### Barcodes
#### Single Cell Multiome Gel Beads A (PN- 2000261)
![image](resources/10x_multinome_beads.png)
* The **10x Barcode** on the ATAC and GEX primers on the same Gel Bead are **NOT identical**. 
Each Gel Bead has a unique pairing of ATAC and GEX barcode. 
* Barcode translation after read processing corresponds ATAC barcodes to GEX barcodes.
* The cell barcode distinguishes between cells, and the UMI (Unique Molecular Identifier) distinguishes between molecules within a cell.

#### Chromium Single Cell Multiome ATAC Library
![image](resources/10x_cellranger_arc_ATAC_barcode.png)
![image](resources/10x_cellranger_arc_ATAC_barcode_seq.png)

#### Chromium Single Cell Multiome GEM Library
![image](resources/10x_cellranger_arc_GEX_barcode.png)
![image](resources/10x_cellranger_arc_GEX_barcode_seq.png)
### Raw fastq
#### ATAC libraries
Consist of standard Illumina® paired-end dsDNA. Sequencing the libraries produces produces a standard Illumina® BCL data 
output folder that includes fastq files: \
`[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz` \
where `[Read Type]` is one of:
* I1: Dual index i7 read (optional) – includes the 8bp Sample Index
* R1: Read 1N
* R2: Dual index i5 read (optional) – includes the 16bp 10x Barcode
* R3: Read 2N

Read 1N and Read 2N are paired-end reads used for sequencing the DNA insert.

#### GEX libraries
Consist of standard Illumina® paired-end cDNA. Sequencing the libraries produces produces a standard Illumina® BCL data 
output folder that includes fastq files: \
`[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz` \
where `[Read Type]` is one of:
* I1: Dual index i7 read (optional) - includes the 10bp sample index sequence i7.
* I2: Dual index i5 read (optional) - includes the 10bp sample index sequence i5.
* R1: TruSeq Read 1 - used to sequence the 16 bp 10x Barcodes + 12 bp UMI
* R2: TruSeq Read 2 - used to sequence the insert

### BAM
#### [ATAC BAM](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/bam-atac)
The `cellranger-arc count` outputs `atac_possorted_bam.bam`, 
a position-sorted and indexed BAM file for the Chromatin Accessibility library.
Tags:
* `CB`	Z   - Chromium **cellular barcode** sequence that is error-corrected, 
confirmed against a list of known-good barcode sequences and translated.(e.g. AGAATGGTCTGCAT-1)
    * Barcode translation - the _in silico_ translation of the error-corrected ATAC barcode
    to its corresponding paired GEX barcode.
* `CR`	Z	- Chromium cellular barcode sequence as reported by the sequencer.
* `BC`	Z	- Sample index read.
* `TR`	Z	- Adapter sequence trimmed off the end of the read.
* `GP`	i	- Genome position. _Note: this is an auxiliary tag used for the purpose of 
duplicate marking and is not intended for downstream use. We intend to deprecate this tag in subsequent versions._

#### [GEX BAM](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/bam-gex)
The `cellranger-arc` outputs a position-sorted and indexed BAM file of read alignments 
to the genome and transcriptome. 
Reads aligned to the transcriptome across **exon junctions** in the genome tend to have a large gap 
in their CIGAR string i.e. 35M225N64M. Each read in this BAM file has a 10x Chromium cellular 
(associated with a 10x gel bead) barcode and molecular barcode information attached. 
Cell Ranger ARC modifies MAPQ values

10x Chromium barcode (associated with a 10x gel bead) and molecular barcode information 
for each read is stored as TAG fields:

* `CB`	Z	- Chromium **cellular barcode** sequence that is error-corrected and confirmed against 
a list of known-good barcode sequences.
* `CR`	Z	- Chromium cellular barcode sequence as reported by the sequencer.
* `UB`	Z	- Chromium **molecular barcode** sequence that is error-corrected among other 
molecular barcodes with the same cellular barcode and gene alignment.
* `UR`	Z	- Chromium molecular barcode sequence as reported by the sequencer.

The following tags will also be present on reads that mapped to the genome and overlapped 
either an exon or an intron by at least one base pair (default mode). 
A read may align to multiple transcripts and genes, but **it is only considered confidently 
mapped to the transcriptome if it maps to a single gene**.

* `RE`	A	Single character indicating the region type of this alignment 
(`E` = exonic, `N` = intronic, `I` = intergenic).
* `TX` 	Z	Semicolon-separated list present in reads aligned to the same strand as the 
transcripts (or genes) compatible with this alignment. In case of a transcriptomic alignment 
overlapping an exonic region the format of each entry is `[transcript_id]`,`[strand][pos]`,`[cigar]`;
 where `transcript_id` is specifed by the reference GTF, strand is either `+` or `-`, `pos` is the 
 alignment offset in transcript coordinates, and `cigar` is the CIGAR string in t
 ranscript coordinates. In case of a genomic alignment overlapping an intronic region the 
 format of each entry is `[gene_id],[strand]`; where `gene_id` is specifed by the reference 
 GTF and strand is either `+` or `-`.
* `AN`	Z	Same as the TX tag, but for reads that are aligned to the antisense strand 
of annotated transcripts (or genes).
* `GX` 	Z	Semicolon-separated list of gene IDs that are compatible with this alignment. 
Gene IDs are specified with the `gene_id` key in the reference GTF attribute column.
* `GN` 	Z	Semicolon-separated list of gene names that are compatible with this alignment. 
Gene names are specified with `gene_name` key in the reference GTF attribute column.
* `MM`	i	Set to 1 if the genome-aligner (STAR) originally gave a **MAPQ < 255** 
(it multi-mapped to the genome) and Cell Ranger changed it to 255 because the read overlapped 
exactly one gene.
* `pa`	i	The **number of poly-A nucleotides trimmed** from the 3' end of read 2. Up to 10% 
mismatches are permitted.
* `ts` 	i	The **number of template switch oligo (TSO) nucleotides trimmed** from 
the 5' end of read 2. Up to 3 mismatches are permitted. The 30-bp TSO sequence 
is `AAGCAGTGGTATCAACGCAGAGTACATGGG`.

#### Transcriptome alignment
<img src="https://support.10xgenomics.com/img/cellranger-arc/gex-processing-introns.png" width=700x>

### [GEX molecule info](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/molecule_info)
An HDF5 file containing per-molecule information for all molecules that contain a valid barcode and valid UMI and were assigned with high confidence to a gene.

HDF5 file hierarchy
```
(root)
├─ barcode_idx
├─ barcode_info [HDF5 group]
│   ├─ genomes
│   └─ pass_filter
├─ barcodes
├─ count
├─ feature_idx
├─ features [HDF5 group]
│   ├─ _all_tag_keys
│   ├─ feature_type
│   ├─ genome
│   ├─ id
│   └─ name
├─ gem_group
├─ library_idx
├─ library_info
├─ metrics_json
├─ umi
└─ umi_type
```
### Per-barcode metrics (csv)
Each row represents a every observed barcode. The columns contain the paired ATAC and Gene Expression barcode sequences, ATAC and Gex QC metrics for that barcode, as well as whether this barcode was identified as a cell-associated partition by the pipeline.
```{bash}
$ head -4 per_barcode_metrics.csv
barcode,gex_barcode,atac_barcode,is_cell,excluded_reason,gex_raw_reads,gex_mapped_reads,gex_conf_intergenic_reads,gex_conf_exonic_reads,gex_conf_intronic_reads,gex_conf_exonic_unique_reads,gex_conf_exonic_antisense_reads,gex_conf_exonic_dup_reads,gex_exonic_umis,gex_conf_intronic_unique_reads,gex_conf_intronic_antisense_reads,gex_conf_intronic_dup_reads,gex_intronic_umis,gex_conf_txomic_unique_reads,gex_umis_count,gex_genes_count,atac_raw_reads,atac_unmapped_reads,atac_lowmapq,atac_dup_reads,atac_chimeric_reads,atac_mitochondrial_reads,atac_fragments,atac_TSS_fragments,atac_peak_region_fragments,atac_peak_region_cutsites
AAACAGCCAAACAACA-1,AAACAGCCAAACAACA-1,ACAGCGGGTGTGTTAC-1,0,0,11,10,1,7,2,6,1,4,2,0,2,0,0,6,2,2,9,0,2,1,0,0,6,4,6,12
AAACAGCCAAACATAG-1,AAACAGCCAAACATAG-1,ACAGCGGGTTGTTCTT-1,0,2,7,7,0,6,1,6,0,5,1,0,1,0,0,6,1,1,3,0,2,0,0,0,1,0,0,0
AAACAGCCAAATATCC-1,AAACAGCCAAATATCC-1,ACAGCGGGTTGTGACT-1,1,0,58263,56138,2528,33577,17962,31480,1827,27942,3431,11276,6479,9899,1316,30809,4747,2272,16029,141,1332,6440,157,22,7937,4822,7074,13986
```

### [TSO (template switch oligo)](https://kb.10xgenomics.com/hc/en-us/articles/360001493051-What-is-a-template-switch-oligo-TSO-)
* The TSO is hybridizes to untemplated C nucleotides added by the reverse transcriptase during 
**reverse transcription**. The TSO **adds a common 5' sequence to full length cDNA** that is used for 
downstream cDNA amplification.
* TSO sequence (30bp): `AAGCAGTGGTATCAACGCAGAGTACATGGG`
* The TSO is used differently between the Single Cell 3’ and Single Cell 5' assay. 
    * In the 3' assay, the polyd(T) is part of the gel bead oligo (+10x Barcode, UMI,  partial Illumina Read 1), with the TSO supplied in the RT Primer.
    * In the 5' assay, the polyd(T) is supplied in the RT Primer, and the TSO is part of the gel bead oligo.
* Tags `ts:i` and `pa:i` in the output GEX BAM files indicate the number of TSO nucleotides 
trimmed from the 5' end of read 2 and the number of poly-A nucleotides trimmed from the 3' end. 
The trimmed bases are present in the sequence of the BAM record and are soft clipped in the CIGAR string.

## [Cell calling](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/algorithms/cell-calling)
Identify the barcodes of cell population from the non-cell background using **paired information**:
* ATAC peak-bc matrix: count of transposition events in peaks
* GEX gene-bc matrix: gene expression UMIs 

Steps:

1. Barcodes filtering
    1. Mask ATAC "low targeting" barcodes: whose **fraction of fragments overlapping called peaks** is **lower than the fraction of genome in peaks** (for this calculation, peaks are padded by 2kbp on both sides to account for fragment length)).

    2. Mask ATAC gel bead doublets: putative gel bead doublets where a partition contains **one cell and two barcoded gel beads**.

    3. Greater than 1 count observed in each library: the minimum threshold is set at >1 count in each library.

2. Cell calling
    1. **De-duplication**
        
        Plots each barcode into a 2D space defined by their **ATAC and GEX counts**. 
        Barcodes with **identical coordinates** are collapsed into a single 
        measurement to generate a more uniform density of points across the 2D space. _(is it reasonable?)_
        
        This de-emphasizes over-represented low-count barcodes and allows 
        suppression of noise without using thresholds or making assumptions about the count 
        distribution profiles.
        
    2. **Ordmag-derived initial grouping** (initial labeling: om_cell/om_non-cell)
        
        Filtering using thresholds derived from "ordmag". Ordmag is a published algorithm that finds 
        **a threshold that is 10-fold less than the maximum value** after removing outliers. 
        A threshold is defined independently for each dimension and barcodes above both ATAC and GEX 
        thresholds are labeled as cells, with the remainder labeled as non-cells.
        
    3. **K-means boundary refinement** (refined labelling: km_cell/km_non-cell)
        
        Using K-means to **refine the boundaries** of these initial set of cells. The K-means is initialized using 
        centroids calculated from the ordmag-defined cell and non-cell groups, setting K=2.
        
    4. **Map classification to de-duplicated barcodes** 
        
        - Step 3 K-means assignment - used to classify the full set of non-excluded barcodes. 
        - Step 1 masked barcodes are assigned to the cell/non-cell based on the K-means classification of their counterpart with identical ATAC and GEX counts.

   **Force cell**: override override the default pipeline Joint Cell Calling algorithm, when additional parameters of `cellranger-arc count` are provided: `--min-atac-count=N`(min number of transposition events in peaks for a cell barcode) and `--min-gex-count=N` (min number of UMI counts 
for a cell barcode).

   <img src="https://support.10xgenomics.com/img/cellranger-arc/joint-cell-calling.png" height="250x"><img src="https://support.10xgenomics.com/img/cellranger-arc/force-cells.png" height="250x">
    

## Downstream Analysis (based on peak-cell count matrix and gene-cell count matrix)
### Reviews
* https://www.nature.com/articles/s41576-019-0093-7
* https://www.nature.com/articles/s41592-019-0692-4#Sec4
* https://academic.oup.com/bib/article/22/1/20/5828125

<img src="https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fs41576-019-0093-7/MediaObjects/41576_2019_93_Fig3_HTML.png?as=webp" width="600x">
* CCA, canonical correlation analysis (identify a set of variables that are maximally correlated between two data sets); 
* MOFA, multi-omics factor analysis (identifying a set of factors that explain variance across multiple data modalities and used their method to jointly analyse bulk genomic);
* NMF, non-negative matrix factorization.

### [10x](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/overview/welcome) 
#### software: `cellranger-arc`
* Dimension reduction
    * PCA for GEX and LSA for ATAC
    * t-SNE
    * UMAP
* Cluster-based Analyses
    * Clustering
    * Differential Enrichment Analysis
* Feature Linkage
    * Transcription Factor Analysis
    * Peak-motif Occurence Mappings
    
### [Seurat](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1)
Use WNN to represent a weighted combination of different modalities. Steps:
* Independent preprocessing and dimensional reduction of each modality individually
* Learning cell-specific modality ‘weights’, and constructing a WNN graph that integrates the modalities
* Downstream analysis (i.e. visualization, clustering, etc.) of the WNN graph


1. pre-processing and dimension reduction on RNA and ATAC:
    ```{r}
    # RNA analysis
    DefaultAssay(pbmc) <- "RNA"
    pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% 
    RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
     ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(pbmc) <- "ATAC"
    pbmc <- RunTFIDF(pbmc)
    pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
    pbmc <- RunSVD(pbmc)
    pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, 
    reduction.name = "umap.atac", reduction.key = "atacUMAP_")```
2. Calculate a WNN graph  -->  UMAP and clustering --> annotate the clusters
3. Visualize clustering based on gene expression, ATAC-seq, or WNN analysis. 
The differences are more subtle than in the previous analysis (**the weights are more evenly split** 
than in our CITE-seq example), but WNN provides the clearest separation of cell states.

    E.g. the ATAC-seq data assists in the separation of CD4 and CD8 T cell states. 
    This is due to the presence of multiple loci that exhibit differential accessibility between 
    different T cell subtypes.
    <img src="https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis_files/figure-html/UMAPs-1.png" width="600x">
4. Examine the accessible regions of each cell to determine enriched motifs. (e.g. ChromVar)
5. Explore the multimodal dataset to identify key regulators (e.g. TF modulators) of each cell state. 

### [ArchR](https://www.archrproject.com/bookdown/co-accessibility-with-archr.html) 
To circumvent substantial noise in correlated analysis caused by sparse single-cell data, 
ArchR adopts an approach introduced by Cicero to create **low-overlapping aggregates 
of single cells**, in which filtering aggregates with greater than 80% 
overlap with any other aggregate in order to reduce bias.

Integrates scATAC-seq and scRNA-seq for “peak-to-gene links” (links the center of the peak to 
the single-base TSS of the gene.) and the prediction of enhancer activity through peak-to-gene linkage analysis.

<img src="https://www.archrproject.com/bookdown/images/HemeWalkthrough/PNG/Plot-Tracks-Marker-Genes-with-Peak2GeneLinks_5.png" width="600x">
<img src="resources/archR_crop_Plot-Tracks-Marker-Genes-with-CoAccessibility_5_cluster2_cd14.png" width="600x">


## Workflow
* snakemake “all_rna”  integrates with rna-seq
   * `snakemake all_rna …`
* Input of gene expression
   10x cellranger-arc output of the raw barcode-feature matrix (we do not want to alter the 10x pre-processing steps for RNA-seq)
* Output
   2 count matrices: cell-peak count matrix; cell-gene expression matrix
* Downstream analysis:
   * co-clustering
   * linkage plot
* Barcodes are selected based on avocato ATAC-seq processing results
* Dataset: 10x multinome pbmc3k granunocytes sorted
   Build: GRCh38
   
* QC on RNA-seq?



