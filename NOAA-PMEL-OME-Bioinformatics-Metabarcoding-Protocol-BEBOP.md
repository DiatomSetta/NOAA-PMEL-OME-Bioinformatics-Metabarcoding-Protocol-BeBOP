---
# MIOP terms
methodology category: Omics Analysis
project: NOAA Pacific Marine Environmental Laboratory Ocean Molecular Ecology Group protocols
purpose: 'taxonomic diversity assessment by targeted gene survey [OBI:0001960]'
analyses: 'bioinformatics analysis [MIOP:0000004] | amplicon sequencing assay [OBI:0002767]'
geographic location: 
   options: 'North East Pacific Ocean [GAZ:00013765] | Bering Sea [GAZ:00008990] | Arctic Ocean [GAZ:00000323]'
broad-scale environmental context: 'marine biome [ENVO:00000447] | oceanic epipelagic zone biome [ENVO:01000035]'
local environmental context: 
   options: 'marine photic zone [ENVO:00000209] | marine aphotic zone [ENVO:00000210] | marine benthic biome [ENVO:01000024]'
environmental medium: 'ocean water [ENVO:00002149] | sea water [ENVO:00002149]'
target: 'Bacteria-16S-V4V5-Parada [ODE]'
creator: Samantha Setta, Sean McAllister, Zachary Gold
materials required: high-performance computing resources
skills required: basic bash, R
time required: variable
personnel required: 1
language: en
issued: TBD
audience: scientists
publisher: 'NOAA Pacific Marine Environmental Laboratory Ocean Molecular Ecology Group | University of Washington Cooperative Institute for Climate, Ocean, & Ecosystem Studies'
hasVersion: 1
license: 'CC0 1.0 Universal'
maturity level: demonstrated

# FAIRe terms
sop_bioinformatics: TBD_this_DOI (paste link when published)
checkls_ver: TBD_bioinformatics_template
assay_name: 'ssu16sv4v5_parada | ssu16sv4v5_parada_OSUmod'
pcr_primer_forward: GTGYCAGCMGCCGCGGTAA
pcr_primer_reverse: CCGYCAATTYMTTTRAGTTT
trim_method: 'Cutadapt, primer trimming | DADA2, filterAndTrim (quality and length trimming)'
trim_param:
   - Cutadapt:
      default: '-a "{primerF};required...{revcomp_primerR};optional", -A "{primerR};required...{revcomp_primerF};optional", --discard-untrimmed, -m 1'
      source_file: REVAMP_config
      source_term: 'primerF | revcomp_primerR | primerR | revcomp_primerF'
   - DADA2:
      default: 'trunQ = {dada_trunQ}, trimRight = {dada_trimRight}, trimLeft = {dada_trimLeft}'
      source_file: REVAMP_config
      source_term: 'dada_trunQ | dada_trimRight | dada_trimLeft'
merge_tool: 'DADA2, mergePairs'
merge_min_overlap: 20
min_len_cutoff:
   default: 100
   source_file: REVAMP_config
   source_term: dada_minlength
min_len_tool: DADA2
error_rate_tool: DADA2
error_rate_type: expected error rate
error_rate_cutoff:
   default: 2
   source_file: REVAMP_config
   source_term: 'dada_maxEE1 | dada_maxEE2'
chimera_check_method: 'DADA2, removeBimeraDenovo, consensus'
chimera_check_param: not applicable
otu_clust_tool: 'DADA2, pool="pseudo"'
otu_clust_cutoff: 100
min_reads_cutoff: 
   default: 2
   options: 'raw = 1 | quality-filtered = 1 | final-filtered >= 2'
   source_file: decontam_workflow_config
   source_term: TBD (n_ton_removal)
min_reads_cutoff_unit: reads
min_reads_tool: 'DADA2 | decontam_workflow'
otu_db:
   - REVAMP:
      default: NCBI GenBank nt database, downloaded {nt_database_version}'
      source_file: REVAMP_config
      source_term: nt_database_version
   - SILVAngs:
      default: 'non-redundant SILVA SSU ref dataset, release {silva_db_release}'
      source_file: REVAMP_config
      source_term: silva_db_release
   - scikit-learn-silva:
      default: 'TBD'
      source_file: TBD
      source_term: TBD
   - Anacapa:
      default: 'rCRUX db {TBD}'
      source_file: TBD
      source_term: TBD
otu_db_custom: not applicable
tax_assign_cat:
   - REVAMP:
      default: sequence similarity
   - SILVAngs:
      default: sequence similarity
   - scikit-learn-silva:
      default: probabilistic
   - Anacapa:
      default: probabilistic
otu_seq_comp_appr:
   - REVAMP:
      default: blastn >2.14.1+
   - SILVAngs:
      default: blastn 2.11.0+
   - scikit-learn-silva:
      default: k-mer based
   - Anacapa:
      default: Bowtie 2
tax_class_id_cutoff:
   - REVAMP:
      default: 60
      options: 'species = 97 | genus = 95 | family = 90 | order = 80 | class = 70 | phylum = 60'
      source_file: REVAMP_config
      source_term: taxonomyConfidenceCutoffs
   - SILVAngs:
      default: 86
   - scikit-learn-silva:
      default: not applicable
   - Anacapa:
      default: 40
      options: 'species = 95 | genus = 90 | family = 80 | order = 70 | class = 60 | phylum = 50 | any = 40'
tax_class_query_cutoff:
   - REVAMP:
      default: 90
      source_file: REVAMP_config
      source_term: blastQueryCovCutoff
   - SILVAngs:
      default: 86
   - scikit-learn-silva:
      default: not applicable
   - Anacapa:
      default: 80
tax_class_collapse:
   - REVAMP:
      default: 'Taxonomic levels were dropped to the lowest common ancestor (LCA), and were further dropped depending on % identity thresholds ({taxonomyConfidenceCutoffs}, see REVAMP readme).'
      source_file: REVAMP_config
      source_term: taxonomyConfidenceCutoffs
   - SILVAngs:
      default: 'No taxonomic levels are dropped from best blast hit matches due to high quality comprehensive nature of the reference dataset and taxonomy.'
   - scikit-learn-silva:
      default: TBD check - 'Taxonomic levels were dropped to the lowest common ancestor (LCA), and were further dropped depending on % confidence thresholds.'
   - Anacapa:
      default: 'Taxonomic levels were dropped to the lowest common ancestor (LCA), and were further dropped depending on % confidence thresholds.'
tax_class_other: not applicable
screen_geograph_method: not applicable
screen_contam_0_1:
   - raw:
      default: 0
   - quality-filtered
      default: 1
   - final-filtered
      default: 1
screen_contam_method: 'Applicable to GBIF/OBIS submission. Not applicable to NODE submission. See below'




screen_nontarget_method: 'Applicable to GBIF/OBIS submission. Not applicable to NODE submission. See below'
screen_other: 'Applicable to GBIF/OBIS submission. Not applicable to NODE submission. See below'
otu_raw_description: 'No filtering outside of DADA2 default ASV denoising'
otu_final_description: this_DOI (link to decontamination screening section)
bioinfo_method_additional: this_DOI (paste link when published)

# NOAA PMEL Ocean Molecular Ecology terms
asv_method
discard_untrimmed: 1
qiime2_version: not applicable
tourmaline_asv_method: not applicable
dada2_trunc_len_f: not applicable
dada2pe_trunc_len_r: not applicable
dada2_trim_left_f: {dada_trimLeft}
dada2pe_trim_left_r: {dada_trimLeft}
dada2_trim_right_f: {dada_trimRight}
dada2pe_trim_right_r: {dada_trimRight}
dada2_max_ee_f: {dada_maxEE1}
dada2pe_max_ee_r: {dada_maxEE2}
dada2_trunc_q: {dada_trunQ}
dada2_pooling_method: pseudo
dada2_chimera_method: 'removeBimeraDenovo, consensus'
dada2_min_fold_parent_over_abundance: not applicable
dada2_n_reads_learn: 2016000000
repseqs_min_abundance: not applicable
repseqs_min_length: not applicable
repseqs_max_length: not applicable
repseqs_min_prevalence: not applicable
min_consensus: not applicable
blca_confidence: not applicable

---

# NOAA PMEL OME Bioinformatics Metabarcoding Protocol

## PROTOCOL INFORMATION

### Minimum Information about an Omics Protocol (MIOP)

- MIOP terms are listed in the YAML frontmatter of this page.
- See <https://github.com/BeBOP-OBON/miop/blob/main/model/schema/terms.yaml> for list and definitions.

### Authors

| PREPARED BY | AFFILIATION | ORCID | DATE |
| ------------- | ------------- | ------------- | ------------- |
| Sean McAllister | Ocean Molecular Ecology, NOAA PMEL & UW CICOES | <https://orcid.org/0000-0001-6654-3495> | 2025-05-08 |
| Samantha Setta | Ocean Molecular Ecology, NOAA PMEL & UW CICOES | <https://orcid.org/0000-0001-9075-7573> | 2025-04-30 |
| Shannon Brown | Ocean Molecular Ecology, NOAA PMEL & UW CICOES | <https://orcid.org/0000-0001-9808-2638> | 2025-04-30 |
| Zachary Gold | Ocean Molecular Ecology, NOAA PMEL | <https://orcid.org/0000-0003-0490-7630> | 2025-04-30 |

- All authors known to have contributed to the preparation of this protocol should be listed, including those who filled in the template.
- Visit <https://orcid.org/> to register for an ORCID.

### Protocol Revision Record

| VERSION | RELEASE DATE | DESCRIPTION OF REVISIONS |
| ------------- | ------------- | ------------- |
| 1.0.0 | 2025-MM-DD | Initial release |

- Version numbers start at 1.0.0 when the protocol is first completed and will increase when changes that impact the outcome of the procedure are made (patches: 1.0.1; minor changes: 1.1.0; major changes: 2.0.0).
- Release date is the date when a given protocol version was finalised.
- Description of revisions includes a brief description of what was changed relative to the previous version.

---

NEED TO EDIT HEADINGS BELOW...

### RELATED PROTOCOLS IN YOUR FOLDER

This is a list of other protocols deposited in your folder which should be known to users of this protocol. For example, if you create a derivative or altered protocol, you would link to the original protocol in the section below. Please include the link to each related protocol. Also include the version number of that protocol when you linked to it.

| PROTOCOL NAME AND LINK  | VERSION The version of the protocol you linked to | RELEASE DATE This is the date corresponding to the version listed to the left |
| ------------- | ------------- | ------------- |
| Content Cell  | Content Cell  | yyyy-mm-dd  |
| Content Cell  | Content Cell  | yyyy-mm-dd  |

### RELATED EXTERNAL PROTOCOLS

This is a list of other protocols that are not in your folder which should be known to users of this protocol. These include, e.g., kit manuals. Please upload all relevant external protocols to Appendix A and link to them here.

| EXTERNAL PROTOCOL NAME AND LINK  | ISSUER / AUTHOR Please note who authored the protocol (this may also be a company name) | ACCESS DATE This is the date you downloaded or scanned the protocol and uploaded it. |
| ------------- | ------------- | ------------- |
| Content Cell  | Content Cell  | yyyy-mm-dd  |
| Content Cell  | Content Cell  | yyyy-mm-dd  |

### Acronyms and Abbreviations

| ACRONYM / ABBREVIATION | DEFINITION |
| ------------- | ------------- |
| 12S | 12S ribosomal nucleic acid sequencing gene region |
| 16Sv4 | 16S ribosomal nucleic acid V4 gene region |
| 18Sv4 | 18S ribosomal nucleic acid V4 gene region |
| 18Sv9 | 18S ribosomal nucleic acid V9 gene region |
| ASV | Amplicon Sequencing Variant |
| CICOES | Cooperative Institute for Climate, Ocean, and Ecosystem Studies |
| COI | Cytochrome c oxidase subunit I gene region |
| ITS1 | Internal Transcribed Spacer 1 region |
| NOAA | National Oceanographic and Atmospheric Administration |
| PMEL | Pacific Marine Environmental Laboratory |
| UW | University of Washington |

### GLOSSARY

| SPECIALISED TERM | DEFINITION |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |

## BACKGROUND

This document describes the required protocol to conduct insert name of the method/protocol.

### Summary

Insert a short description of the background for the method/protocol (e.g. why and for which purpose do you perform water sampling).
Please provide a brief summary of your method including, as appropriate, a brief description of what techniques your best practice is about, which ocean environments or regions it targets, the primary sensors covered, what type of data/measurements/observing platform it covers, limits to its applicability.

### Method description and rationale

Insert a short description of the functioning principal of the methodology used in the protocol (i.e. how does the method work?). Please note that this is different from the step-by-step description of the protocol procedure.
Insert a short statement explaining why the specific methodology used in the protocol has been selected (e.g. it is highly reproducible, highly accurate, procedures are easy to execute etc….).

Data is run through this standard operating procedure at the sequencing run level, including all samples from the run, regardless of source project. The purpose of this is to maintain reproducibility and operate DADA2 as intended, allowing for proper error learning.

## STANDARD OPERATING PROCEDURE

### Raw Data Download and QA/QC

Raw reads (fastq.gz) are provided by the sequencing center demultiplexed by sample and marker (see sequencing center BeBops listed in RELATED PROTOCOLS). On downloading reads to the local compute infrastructure, fidelity of the downloaded file is checked via a `md5sum` check. Md5s are either supplied by the sequencing center with the original download (fastq.gz.md5) or are provided separately.

To check the quality of the sequencing run, [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is run on all fastq files (`~/path/to/FastQC/fastqc -o fastqc -t 10 *gz`). In addition for checking the run for a number of quality statistics, the "Per base sequence quality" is used to determine the best cutoff for the `dada_trimRight` (number of bases to trim from right or end of sequence) and `dada_trimLeft` (number of bases to trim from the left or start of sequence) parameters for the REVAMP configuration file.

### Rapid Exploration and Visualization through an Automated Metabarcoding Pipeline (REVAMP)

[REVAMP](https://github.com/McAllister-NOAA/REVAMP) ([McAllister et al., 2023](https://doi.org/10.5670/oceanog.2023.231)) is a published metabarcoding pipeline that integrates several tools for read processing, amplicon sequence variant (ASV) assignment, taxonomic assignment, and basic data visualization. The front end of REVAMP, including read trimming and ASV assignment, is used for all markers to produce a consistent set of ASVs to work from for downstream analysis. Taxonomic assignment in REVAMP is also used for downstream analysis on all markers, though some have additional alternative taxonomic classification methods (detailed below).














### Taxonomic Classification

Taxonomic classification occurs after ASV assignment and varies by metabarcoding marker region. 

The 18Sv4 rRNA, 18Sv9 rRNA, 16Sv4 rRNA, and ITS1 regions are classified using [Qiime2 (v.2024.10)](https://qiime2.org/) (Bolyen et al., 2019) feature classifier's [naive bayesian classifier scikit-learn](https://scikit-learn.org/stable/modules/naive_bayes.html). Sci-kit learn classifiers are trained using the [PR2 database(v5.1.0)](https://pr2-database.org/) (Guillou et al., 2012) for 18Sv4 rRNA and 18Sv9 rRNA, [silva (v138.99)](https://docs.qiime2.org/2022.11/data-resources/) for 16Sv4 rRNA (Quast et al., 2013), and a [custom curated database](https://zenodo.org/records/15351664) for the ITS1 region. Each database is curated, extracted to the region of interest, and used to train a naive bayesian taxonomic classifier by:

1) Extracting the region of interest from all sequences in the database with primers used for amplification (see [NOAA PMEL OME Github](https://github.com/NOAA-PMEL/Ocean-Molecular-Ecology)).
```
qiime feature-classifier extract-reads
# --p-min-length and --p-max-length vary by region and are included in the config file.
```
2) Adding sequences that contain the region of interest but do not contain the primer region by alignment to the extracted sequence database. This step is done twice for some metabarcoding regions until ~50-70% of query sequences are retained following [recommendations](https://forum.qiime2.org/t/using-rescripts-extract-seq-segments-to-extract-reference-sequences-without-pcr-primer-pairs/23618) for qiime2 scikit-learn classifiers (see config file).
```
qiime rescript extract-seq-segments
# --p-perc-identity and --p-min-seq-len vary by region and are included in the config file.
```
3) Taxa of sequences no longer in the database are filtered out.
```
qiime rescript filter-taxa
```
4) Sequences are dereplicated to remove any duplicates.
```
qiime rescript dereplicate
# --p-mode uniq
```
5) The taxonomic classifier for each region of interest is trained using qiime2's feature-classifier.
```
qiime feature-classifier fit-classifier-naive-bayes
```
6) The trained classifiers are used to classify ASVs using qiime2's feature classifier.
```
qiime tools export
```

### Decontamination steps

Decontamination occurs after assigning ASVs with revamp and Dada2, to remove ASVs with too few reads or obvious contaminants. The [decontam package](https://doi.org/10.1186/s40168-018-0605-2) (Davis et al., 2018) is used to filter out ASVs, with the following steps:

**1) Estimation of Tag-jumping or Index Hopping.**

Subtract the proportion of reads that jumped into control samples from each environmental sample. Determine which ASVs came from controls vs environmental samples, create a vector of ASVs in positive controls, calculate what proportion of the reads found in the positive controls are found in environment samples, and subtract the composition of the positive controls from the environment samples.
   
**2) Discarding samples with low number of reads.**

Discard samples with less than 10,000 reads, to filter out samples with low sequencing depth.

**3) Clearance of Negative Control Contaminants.**

Identify contaminant ASVs that occurred in field, extraction, and PCR negative controls using the [microDecon package](https://github.com/donaldtmcknight/microDecon). Compare the prevalence of ASVs in blanks to those in samples and remove contaminant ASVs.
   
**4) Remove obvious contaminants (e.g. human, rat, cat, dog).**

Remove any ASVs identified as human (*Homo sapiens*), dog (*Canis familiaris*), cat (*Felis catus*), or other common obvious contaminants.

**5) Remove samples that have one ASV making up 99% of total sequences.**

Remove samples dominated (>99% abundance) by one ASV, suggesting an error with sample processing or sequencing.

# PERSONNEL REQUIRED

Insert the number of technicians, data managers, and scientists required for the good execution of the procedure

## Safety

Identify hazards associated with the procedure and specify protective equipment and safety training required to safely execute the procedure

## Training requirements

Specify technical training required for the good execution of the procedure.

## Time needed to execute the procedure

Specify how much time is necessary to execute the procedure.

# EQUIPMENT, SOFTWARE & PACKAGES


| NAME | VERSION OR MODEL | MANUFACTURER OR CREATOR | REMARKS |
| ------------- | ------------- | ------------- | ------------- |
| Equipment |
| e.g. Laptop | Content Cell | Content Cell | e.g. needs at least 16 GB of RAM |
| Content Cell | Content Cell | Content Cell | Content Cell |
| Software |
| Content Cell | Content Cell | Content Cell | Content Cell |
| Content Cell | Content Cell | Content Cell | Content Cell |
| Code |
| Please include the links to the code you used for this analysis |
| e.g. link to the released version of a github repository  | Content Cell | Content Cell | Content Cell |
| Content Cell | Content Cell | Content Cell | Content Cell |

# GUIDE TO ARCHIVED METHODOLOGY

The contents of this archive should allow your analysis to be reproduced exactly as you intended it.

This document provides guidance on the contents of each partner's compressed archive of in-silico methods. This document should be part of that same archive, serving as an extended README.

Below, please find guidance on what this archive should include. When describing the contents of the archive, please give precise file names and relative paths to the files.




# Archive content

To reproduce the in-silico analysis, please provide one of the following (in order of decreasing preference)

1. Jupyter, R notebook(s) or equivalents

2. Downloaded archive of (the released version of) your github repository

3. Individual scripts

In each of the above cases, guidance and documentation for all the steps you took to perform the in-silico analysis should be included. In case 1., code and documentation are integrated. In cases 2. and 3., in-line comments may be provided, however, these are not generally sufficient as documentation. In those cases, please provide a step-by-step protocol on how and when to run each script in the Execution Procedure section below.

Please include a script on **data acquisition** (e.g. documentation and code to pull sequences from INSDC, access sequences on an institutional FTP server, download metadata files, check file integrity via md5 checksum). Please add sufficient detail, so that the partners only have to install the software, run this script and will then have all the data needed to perform any analysis described below.

### Code
Here please describe each file containing code, including its purpose, its input, its output. Please provide the names and the relative paths to this documentation.

### Code documentation
Here, please indicate if your documentation is with the code (in a code notebook) or stored separately. In-line comments are not considered documentation. If the documentation is stored separately, please provide the names and the relative paths to this documentation.

### Metadata
Please provide link(s) to the files containing metadata about your sequence data (e.g. environmental data, procedural data). Please see the MIxS compliant metadata guidance.

Auxiliary files
e.g. mapping files, test/dummy files, colour palette


## Execution Procedure

Please fill out this section if you have not already documented it as part of your R, Jupyter, or similar notebook. In this section, please provide a step-by-step guidance on how and when to run each component of your code.

## Quality control

In this section please include the names and paths that can be used to validate that operations were successful. If such checks were done during the execution procedures, please note this here. We recommend identifying such steps with in-line tags (e.g. “#QC”).

## Basic troubleshooting guide

Identify known issues associated with the procedure, if any.
Provide troubleshooting guidelines when available.

# REFERENCES

Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019). [https://doi.org/10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9)

Davis, N.M., Proctor, D.M., Holmes, S.P. Relman, D.A., & Callahan, B.J. (2018). Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6, 226. [https://doi.org/10.1186/s40168-018-0605-2](https://doi.org/10.1186/s40168-018-0605-2)

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., ... & Christen, R. (2012). The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. Nucleic acids research, 41(D1), D597-D604. [https://doi.org/10.1093/nar/gks1160](https://doi.org/10.1093/nar/gks1160)

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Opens external link in new windowNucl. Acids Res. 41 (D1): D590-D596.[https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)

DADA2 Callahan et al., 2016

Cutadapt (Martin, 2011)

REVAMP ([McAllister et al., 2023]

# APPENDIX A: DATASHEETS

Link to any documents such as software guidelines, images, etc that support this protocol. Please include a short note describing the document's relevance.
