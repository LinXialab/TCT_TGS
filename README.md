TE Chimeric Transcripts Identification Pipeline
This pipeline is designed for the identification and analysis of transposable element (TE) chimeric transcripts (TCTs) in tumor cells using long-read RNA sequencing. The workflow consists of the following key steps:

Alignment of Long-Read RNA-Seq Data
Long-read RNA sequencing data are aligned to the GRCh38 human genome reference using minimap2 (v2.2.17) with the parameters -ax splice -secondary=no to ensure accurate alignment of splice junctions.

Transcript Assembly and Annotation
Transcript assembly is performed using FLAIR (v1.7.0), followed by FLAIR-correct to refine splice site alignments based on GENCODE v.42 annotations. FLAIR-collapse is applied to generate high-confidence transcript sets, supported by at least ten reads with a MAPQ > 10.
A custom annotation pipeline identifies transcripts that overlap with RepeatMasker-annotated TEs and GENCODE v.42 transcripts. Candidate TCTs are selected based on the presence of at least five supporting reads that span both TEs and adjacent exons of target genes.

Generation of a Reference Transcriptome Including TCTs
The StringTie merge function is used to generate a comprehensive reference transcriptome, integrating GENCODE v.42 annotations with TE-chimeric events that meet predefined filtering criteria.

Transcript-Level Quantification and Candidate Selection
To quantify transcript expression, StringTie’s quantification function (-e -b) is applied to the merged transcriptome. FPKM values are extracted from Ballgown output files to enable transcript-level expression analysis, providing a detailed view of the contribution of TCTs to overall gene expression.

Open-Reading Frame Prediction and Annotation
Open-reading frames (ORFs) within TCTs are predicted using TransDecoder (v5.5.0). The resulting FASTA file contains potential transcript products, classified based on predicted start codon status. TCTs are categorized as normal, truncated, chimeric normal, frame-shift, or chimeric truncated, depending on their coding potential.

This pipeline facilitates the comprehensive identification and functional annotation of TE-chimeric transcripts, offering valuable insights into their potential role in cancer biology.
