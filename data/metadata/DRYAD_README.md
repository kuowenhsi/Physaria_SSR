# README

## Dataset title
Data supporting microsatellite marker development and genotyping in *Physaria globosa*

## Associated manuscript
These files support the manuscript:

**A disjunct distribution and population fragmentation shape rangewide genetic diversity and structure of the endangered Physaria globosa (Brassicaceae)**  
**Christine E. Edwards, Wen-Hsi Kuo, Clara Landon, Burgund Bassüner, Alex G. Linan, and Matthew A. Albrecht**  
**[Insert journal / preprint / DOI here]**

## Overview
This dataset contains files associated with the development and application of 20 microsatellite markers for *Physaria globosa*. The dataset includes:

1. raw Illumina MiSeq paired-end reads from two *P. globosa* individuals used for microsatellite discovery,
2. De novo assembled consensus sequences used to identify microsatellite loci and design primers, and
3. raw capillary electrophoresis trace files used for microsatellite genotyping.

In this study, 16 previously developed microsatellite markers from *Physaria filiformis* were screened for amplification in *P. globosa*. Eleven amplified successfully, and 10 were polymorphic and retained for genotyping. Additional loci were developed specifically for *P. globosa* using shotgun sequencing of genomic DNA from two individuals, de novo assembly, and microsatellite discovery from assembled contigs. From these, 10 additional polymorphic loci were identified. In total, 20 loci were used to genotype 406 samples.

## File inventory

### 1. `P-glo-1_S8_L001_R1_001.fastq`
Raw Illumina MiSeq forward reads (R1) for a first *Physaria globosa* individual used for microsatellite discovery.

- File type: FASTQ
- Sequencing platform: Illumina MiSeq
- Read type: paired-end
- Read length: 150 bp
- Sample role: microsatellite marker development

### 2. `P-glo-1_S8_L001_R2_001.fastq`
Raw Illumina MiSeq reverse reads (R2) for the same *Physaria globosa* individual as above.

- File type: FASTQ
- Sequencing platform: Illumina MiSeq
- Read type: paired-end
- Read length: 150 bp
- Sample role: microsatellite marker development

### 3. `P-glo-2_S9_L001_R1_001.fastq`
Raw Illumina MiSeq forward reads (R1) for a second *Physaria globosa* individual used for microsatellite discovery.

- File type: FASTQ
- Sequencing platform: Illumina MiSeq
- Read type: paired-end
- Read length: 150 bp
- Sample role: microsatellite marker development

### 4. `P-glo-2_S9_L001_R2_001.fastq`
Raw Illumina MiSeq reverse reads (R2) for the same *Physaria globosa* individual as above.

- File type: FASTQ
- Sequencing platform: Illumina MiSeq
- Read type: paired-end
- Read length: 150 bp
- Sample role: microsatellite marker development

### 5. `physaria_Consensus_Sequences.fasta`
Consensus sequences generated from de novo assembly of Illumina MiSeq reads from the two *Physaria globosa* individuals used for marker development. These contig sequences were used for microsatellite identification and primer design.

- File type: FASTA
- Content: assembled contig consensus sequences
- Downstream use: microsatellite discovery with MSATCOMMANDER and primer design with PRIMER3

### 6. `Microsatellite_raw_trace_fsa.zip`
Compressed archive of raw capillary electrophoresis trace files (`.fsa`) generated during microsatellite genotyping.

- File type: ZIP archive containing ABI `.fsa` files
- Platform: Applied Biosystems 3730xl
- Downstream use: allele scoring in GeneMarker v2.6.2
- Content: raw fragment analysis traces for microsatellite loci

## Methods summary

### Marker screening from previously developed loci
We first tested 16 microsatellite markers previously developed for *Physaria filiformis* for amplification in *P. globosa*. PCR amplifications were performed in 10 μL reactions containing GoTaq Flexi DNA polymerase (Promega), 1× Colorless GoTaq Flexi Buffer, 1.5 mM MgCl2, reverse primer, an M13-tagged forward primer, one of four fluorescently labeled M13 primers (6-FAM, VIC, NED, or PET), and dNTPs. PCR products were first screened by agarose gel electrophoresis. Eleven loci amplified successfully in *P. globosa*, and 10 of those were polymorphic and retained for genotyping.

### Marker development from shotgun sequencing
To develop additional markers specific to *P. globosa*, genomic DNA from two individuals was prepared with Nextera DNA sample prep kits and Nextera index kits (Illumina) and sequenced on an Illumina MiSeq using 2 × 150 bp reads. Reads were trimmed and assembled de novo in Geneious v7.1.9 using the “Medium sensitivity/fast” setting. Microsatellites were identified with MSATCOMMANDER, and primers were designed with PRIMER3. An M13 tag was added to the 5′ end of each forward primer to enable universal fluorescent labeling. From the candidate loci screened, 10 additional polymorphic loci with reliable amplification in *P. globosa* were retained.

### Genotyping and allele scoring
All genotyping was carried out by capillary electrophoresis on an Applied Biosystems 3730xl. Electropherograms were scored in GeneMarker v2.6.2 using automated scoring panels developed for each locus and then checked manually.

## Software and tools
The following software and tools were used in generating or processing these data:

- Geneious v7.1.9
- MSATCOMMANDER
- PRIMER3
- GeneMarker v2.6.2
- Applied Biosystems 3730xl capillary sequencer
- Illumina MiSeq

## File format notes

### FASTQ files
The four `.fastq` files are raw paired-end Illumina reads.  
- `R1` = forward reads  
- `R2` = reverse reads

### FASTA file
The `.fasta` file contains assembled consensus/contig sequences derived from the MiSeq reads and used for microsatellite discovery.

### FSA files
The `.fsa` files are raw fragment analysis trace files from capillary electrophoresis. These files can be opened with GeneMarker and other compatible fragment analysis software.

## Contact
For questions about these data, please contact:

Christine E. Edwards
Christy Edwards <Christine.Edwards@mobot.org>
