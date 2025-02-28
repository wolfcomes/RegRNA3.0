# RegRNA 3.0: A Comprehensive Platform for Regulatory RNA Analysis

RegRNA 3.0 is a powerful meta-server and computational workflow designed for the identification and analysis of regulatory RNA motifs, interactions, and annotations. It integrates **26 specialized algorithms** and automated data retrieval from **28 databases** to provide accurate and customizable RNA analysis. This version introduces advanced features such as RNA 3D structure prediction, non-coding RNA (ncRNA) detection, RNA modification analysis, and interaction prediction, making it a versatile tool for researchers in RNA biology and functional genomics.

**Access RegRNA 3.0 at**: [http://awi.cuhk.edu.cn/~RegRNA](http://awi.cuhk.edu.cn/~RegRNA)

---

## Table of Contents
- [Key Features](#key-features)
- [Integrated Tools and Databases](#integrated-tools-and-databases)
  - [RNA Functional Motif Prediction](#rna-functional-motif-prediction)
  - [RNA Interaction Motif Prediction](#rna-interaction-motif-prediction)
  - [RNA Annotation Modules](#rna-annotation-modules)
- [User Interface and Workflow](#user-interface-and-workflow)
- [Installation and Usage](#installation-and-usage)
- [Citation](#citation)
- [Data Availability](#data-availability)
- [Contact](#contact)

---

## Key Features
- **RNA Functional Motif Prediction**: Identify splice sites, riboswitches, AU-rich elements, G-quadruplexes, and more.
- **RNA Interaction Analysis**: Predict miRNA targets, RNA-protein interactions, RNA-ligand binding, and transcription factor binding sites.
- **Comprehensive RNA Annotation**: Classify RNA families, predict subcellular localization, identify RNA modifications, and predict 3D structures.
- **Enhanced Visualization**: Interactive maps, tables, and structural diagrams for intuitive result interpretation.
- **User-Friendly Interface**: Streamlined workflow with parallel processing and customizable parameters.

---

## Integrated Tools and Databases

### RNA Functional Motif Prediction
RegRNA 3.0 integrates a wide range of tools and databases for identifying functional RNA motifs:
- **Splicing Sites**: GeneSplicer ([29](#references)).
- **Polyadenylation Sites**: polya_svm ([30](#references)).
- **Ribosome Binding Sites**: RBSfinder ([31](#references)).
- **Rho-Independent Terminators**: TransTermHP ([32](#references)).
- **AU-Rich Elements**: ARED-Plus ([5](#references)) and ARESite2 ([6](#references)).
- **Riboswitches**: Ribocentre-switch ([7](#references)), RiboSW ([8](#references)), and Rfam ([9](#references)).
- **Core Promoter Elements**: ElemeNT2023 ([11](#references)).
- **RNA Decay**: RNAdegformer ([34](#references)) and OpenVaccine dataset ([12](#references)).
- **G-Quadruplexes**: QUADRatlas ([13](#references)) and G4 Hunter ([47](#references)).

### RNA Interaction Motif Prediction
RegRNA 3.0 supports the prediction of RNA interactions with other molecules:
- **miRNA Targets**: miRanda ([37](#references)) with miRBase ([15](#references)) and miRTarBase ([16](#references)).
- **ncRNA Hybridization Sites**: BLAST ([35](#references)) and RNAcofold ([40](#references)).
- **Transcription Factor Binding**: TRANSFAC ([18](#references)) and Match ([38](#references)).
- **RNA-Ligand Interactions**: RNALigand ([39](#references)) with PDB ([19](#references)), R-BIND ([20](#references)), and Inforna ([21](#references)).
- **RNA-Binding Proteins**: BRIO ([22](#references)).

### RNA Annotation Modules
RegRNA 3.0 provides comprehensive annotation tools for RNA sequences:
- **RNA Families**: RNAcentral ([23](#references)) and eRNABase ([24](#references)).
- **Blood Exosome RNAs**: exoRBase 2.0 ([25](#references)).
- **Subcellular Localization**: RNALocate 3.0 ([26](#references)).
- **A-to-I RNA Editing**: REDIportal 2.0 ([27](#references)).
- **RNA Modifications**: RMVar 2.0 ([28](#references)) and MODOMICS ([41](#references)).
- **RNA 3D Structures**: RhoFold ([42](#references)).
- **Secondary Structures**: RNAfold ([40](#references)), e-RNA ([43](#references)), and KnotFold ([44](#references)).

---

## User Interface and Workflow
RegRNA 3.0 offers a user-friendly interface with the following key components:
1. **Input Section**: Submit RNA sequences in FASTA format (manual entry or file upload).
2. **Analysis Modules**:
   - **Functional Motifs**: Select from 17 tools to analyze splicing, decay, and regulatory elements.
   - **Interaction Motifs**: Predict miRNA targets, RNA-protein interactions, and more.
   - **Annotation**: Classify RNA families, predict localization, and identify modifications.
3. **Results Visualization**:
   - **Map View**: Interactive motif maps with positional annotations.
   - **Table View**: Detailed tables of predicted motifs and their properties.
   - **Structure Displays**: 2D and 3D RNA structure visualizations (forward/reverse orientation).
4. **Export Options**: Download results in tab-delimited or XML formats.

![RegRNA 3.0 Workflow](media/image2.png)
*Figure 1: RegRNA 3.0 workflow integrating data resources, prediction tools, and visualization modules.*

---

## Installation and Usage
RegRNA 3.0 is a web-based tool and does not require installation. Simply visit the web server and follow these steps:
1. Go to [http://awi.cuhk.edu.cn/~RegRNA](http://awi.cuhk.edu.cn/~RegRNA).
2. Input your RNA sequence in FASTA format or upload a file.
3. Select the desired analysis modules and configure parameters.
4. Submit the job and view the results in the interactive interface.

---

## Citation
If you use RegRNA 3.0 in your research, please cite:  
DOI: [insert DOI upon publication].

---

## Data Availability
- **RegRNA 3.0 Web Server**: [http://awi.cuhk.edu.cn/~RegRNA](http://awi.cuhk.edu.cn/~RegRNA)
- **Supplementary Data**: Available online at NAR.

---

## Contact
For technical support or inquiries, please contact:  
**Email**: [insert correspondence email]  
**Developed by**: [Author Names and Affiliations]

---

## References
For a full list of databases and tools integrated into RegRNA 3.0, refer to the [References](#references) section in the manuscript.