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

![RegRNA 3.0 Workflow](Figures/F3.png)
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

## Contact
For technical support or inquiries, please contact:  
**Email**: [insert correspondence email]  
**Developed by**: [Author Names and Affiliations]

---

## References
## References
1. Stamm, S., Riethoven, J.-J., Le Texier, V., Gopalakrishnan, C., Kumanduri, V., Tang, Y., Barbosa-Morais, N.L. and Thanaraj, T.A. (2006) ASD: a bioinformatics resource on alternative splicing. *Nucleic Acids Research*, **34**, D46-D55.
2. Giulietti, M., Piva, F., D'Antonio, M., D'Onorio De Meo, P., Paoletti, D., Castrignano, T., D'Erchia, A.M., Picardi, E., Zambelli, F. and Principato, G. (2013) SpliceAid-F: a database of human splicing factors and their RNA-binding sites. *Nucleic Acids Research*, **41**, D125-D131.
3. Grillo, G., Turi, A., Licciulli, F., Mignone, F., Liuni, S., Banfi, S., Gennarino, V.A., Horner, D.S., Pavesi, G. and Picardi, E. (2010) UTRdb and UTRsite (RELEASE 2010): a collection of sequences and regulatory motifs of the untranslated regions of eukaryotic mRNAs. *Nucleic Acids Research*, **38**, D75-D80.
4. Lo Giudice, C., Zambelli, F., Chiara, M., Pavesi, G., Tangaro, M.A., Picardi, E. and Pesole, G. (2023) UTRdb 2.0: a comprehensive, expert curated catalog of eukaryotic mRNAs untranslated regions. *Nucleic Acids Research*, **51**, D337-D344.
5. Bakheet, T., Hitti, E. and Khabar, K.S.A. (2018) ARED-Plus: an updated and expanded database of AU-rich element-containing mRNAs and pre-mRNAs. *Nucleic Acids Research*, **46**, D218-D220.
6. Fallmann, J., Sedlyarov, V., Tanzer, A., Kovarik, P. and Hofacker, I.L. (2016) AREsite2: an enhanced database for the comprehensive investigation of AU/GU/U-rich elements. *Nucleic Acids Research*, **44**, D90-D95.
7. Bu, F., Lin, X., Liao, W., Lu, Z., He, Y., Luo, Y., Peng, X., Li, M., Huang, Y. and Chen, X. (2024) Ribocentre-switch: a database of riboswitches. *Nucleic Acids Research*, **52**, D265-D272.
8. Chang, T.-H., Huang, H.-D., Wu, L.-C., Yeh, C.-T., Liu, B.-J. and Horng, J.-T. (2009) Computational identification of riboswitches based on RNA conserved functional sequences and conformations. *RNA*, **15**, 1426-1430.
9. Ontiveros-Palacios, N., Cooke, E., Nawrocki, E.P., Triebel, S., Marz, M., Rivas, E., Griffiths-Jones, S., Petrov, A.I., Bateman, A. and Sweeney, B. (2024) Rfam 15: RNA families database in 2025. *Nucleic Acids Research*, gkae1023.
10. Lambert, A., Fontaine, J.-F., Legendre, M., Leclerc, F., Permal, E., Major, F., Putzer, H., Delfour, O., Michot, B. and Gautheret, D. (2004) The ERPIN server: an interface to profile-based RNA motif identification. *Nucleic Acids Research*, **32**, W160-W165.
11. Adato, O., Sloutskin, A., Komemi, H., Brabb, I., Duttke, S., Bucher, P., Unger, R. and Juven-Gershon, T. (2024) ElemeNT 2023: an enhanced tool for detection and curation of core promoter elements. *Bioinformatics*, **40**, btae110.
12. Wayment-Steele, H.K., Kladwang, W., Watkins, A.M., Kim, D.S., Tunguz, B., Reade, W., Demkin, M., Romano, J., Wellington-Oguri, R. and Nicol, J.J. (2022) Deep learning models for predicting RNA degradation via dual crowdsourcing. *Nature Machine Intelligence*, **4**, 1174-1184.
13. Bourdon, S., Herviou, P., Dumas, L., Destefanis, E., Zen, A., Cammas, A., Millevoi, S. and Dassi, E. (2023) QUADRatlas: the RNA G-quadruplex and RG4-binding proteins database. *Nucleic Acids Research*, **51**, D240-D247.
14. Mituyama, T., Yamada, K., Hattori, E., Okida, H., Ono, Y., Terai, G., Yoshizawa, A., Komori, T. and Asai, K. (2009) The Functional RNA Database 3.0: databases to support mining and annotation of functional RNAs. *Nucleic Acids Research*, **37**, D89-D92.
15. Kozomara, A., Birgaoanu, M. and Griffiths-Jones, S. (2019) miRBase: from microRNA sequences to function. *Nucleic Acids Research*, **47**, D155-D162.
16. Cui, S., Yu, S., Huang, H.-Y., Lin, Y.-C.-D., Huang, Y., Zhang, B., Xiao, J., Zuo, H., Wang, J. and Li, Z. (2024) miRTarBase 2025: updates to the collection of experimentally validated microRNA--target interactions. *Nucleic Acids Research*, gkae1072.
17. Zhao, Y., Li, H., Fang, S., Kang, Y., Wu, W., Hao, Y., Li, Z., Bu, D., Sun, N. and Zhang, M.Q. (2016) NONCODE 2016: an informative and valuable data source of long non-coding RNAs. *Nucleic Acids Research*, **44**, D203-D208.
18. Matys, V., Kel-Margoulis, O.V., Fricke, E., Liebich, I., Land, S., Barre-Dirrie, A., Reuter, I., Chekmenev, D., Krull, M. and Hornischer, K. (2006) TRANSFAC® and its module TRANSCompel®: transcriptional gene regulation in eukaryotes. *Nucleic Acids Research*, **34**, D108-D110.
19. Burley, S.K., Berman, H.M., Kleywegt, G.J., Markley, J.L., Nakamura, H. and Velankar, S. (2017) Protein Data Bank (PDB): the single global macromolecular structure archive. *Protein Crystallography: Methods and Protocols*, 627-641.
20. Donlic, A., Swanson, E.G., Chiu, L.-Y., Wicks, S.L., Juru, A.U., Cai, Z., Kassam, K., Laudeman, C., Sanaba, B.G. and Sugarman, A. (2022) R-BIND 2.0: an updated database of bioactive RNA-targeting small molecules and associated RNA secondary structures. *ACS Chemical Biology*, **17**, 1556-1566.
21. Disney, M.D., Winkelsas, A.M., Velagapudi, S.P., Southern, M., Fallahi, M. and Childs-Disney, J.L. (2016) Inforna 2.0: a platform for the sequence-based design of small molecules targeting structured RNAs. *ACS Chemical Biology*, **11**, 1720-1728.
22. Guarracino, A., Pepe, G., Ballesio, F., Adinolfi, M., Pietrosanto, M., Sangiovanni, E., Vitale, I., Ausiello, G. and Helmer-Citterich, M. (2021) BRIO: a web server for RNA sequence and structure motif scan. *Nucleic Acids Research*, **49**, W67-W71.
23. (2021) RNAcentral 2021: secondary structure integration, improved sequence search and new member databases. *Nucleic Acids Research*, **49**, D212-D220.
24. Song, C., Zhang, G., Mu, X., Feng, C., Zhang, Q., Song, S., Zhang, Y., Yin, M., Zhang, H. and Tang, H. (2024) eRNAbase: a comprehensive database for decoding the regulatory eRNAs in human and mouse. *Nucleic Acids Research*, **52**, D81-D91.
25. Lai, H., Li, Y., Zhang, H., Hu, J., Liao, J., Su, Y., Li, Q., Chen, B., Li, C. and Wang, Z. (2022) exoRBase 2.0: an atlas of mRNA, lncRNA and circRNA in extracellular vesicles from human biofluids. *Nucleic Acids Research*, **50**, D118-D128.
26. Wu, L., Wang, L., Hu, S., Tang, G., Chen, J., Yi, Y., Xie, H., Lin, J., Wang, M. and Wang, D. (2024) RNALocate v3.0: Advancing the Repository of RNA Subcellular Localization with Dynamic Analysis and Prediction. *Nucleic Acids Research*, gkae872.
27. Mansi, L., Tangaro, M.A., Lo Giudice, C., Flati, T., Kopel, E., Schaffer, A.A., Castrignanò, T., Chillemi, G., Pesole, G. and Picardi, E. (2021) REDIportal: millions of novel A-to-I RNA editing events from thousands of RNAseq experiments. *Nucleic Acids Research*, **49**, D1012-D1019.
28. Huang, Y., Zhang, L., Mu, W., Zheng, M., Bao, X., Li, H., Luo, X., Ren, J. and Zuo, Z. (2024) RMVar 2.0: an updated database of functional variants in RNA modifications. *Nucleic Acids Research*, gkae924.
29. Pertea, M., Lin, X. and Salzberg, S.L. (2001) GeneSplicer: a new computational method for splice site prediction. *Nucleic Acids Research*, **29**, 1185-1190.
30. Cheng, Y., Miura, R.M. and Tian, B. (2006) Prediction of mRNA polyadenylation sites by support vector machine. *Bioinformatics*, **22**, 2320-2325.
31. Suzek, B.E., Ermolaeva, M.D., Schreiber, M. and Salzberg, S.L. (2001) A probabilistic method for identifying start codons in bacterial genomes. *Bioinformatics*, **17**, 1123-1130.
32. Kingsford, C.L., Ayanbule, K. and Salzberg, S.L. (2007) Rapid, accurate, computational discovery of Rho-independent transcription terminators illuminates their relationship to DNA uptake. *Genome Biology*, **8**, 1-12.
33. Du, P. and Li, Y. (2008) Prediction of C-to-U RNA editing sites in plant mitochondria using both biochemical and evolutionary information. *Journal of Theoretical Biology*, **253**, 579-586.
34. He, S., Gao, B., Sabnis, R., and Sun, Q. (2023). RNAdegformer: accurate prediction of mRNA degradation at nucleotide resolution with deep learning. *Briefings in bioinformatics*, **24**, bbac581.
35. Altschul, S.F., Gish, W., Miller, W., Myers, E.W., and Lipman, D.J. (1990). Basic local alignment search tool. *Journal of molecular biology*, **215**, 403-410.
36. Grillo, G., Licciulli, F., Liuni, S., Sbisa, E., and Pesole, G. (2003). PatSearch: a program for the detection of patterns and structural motifs in nucleotide sequences. *Nucleic acids research*, **31**, 3608-3612.
37. John, B., Enright, A.J., Aravin, A., Tuschl, T., Sander, C., and Marks, D.S. (2004). Human microRNA targets. *PLoS biology*, **2**, e363.
38. Kel, A.E., Gossling, E., Reuter, I., Cheremushkin, E., Kel-Margoulis, O.V., and Wingender, E. (2003). MATCHTM: a tool for searching transcription factor binding sites in DNA sequences. *Nucleic acids research*, **31**, 3576-3579.
39. Sun, S., Yang, J., and Zhang, Z. (2022). RNALigands: a database and web server for RNA–ligand interactions. *RNA*, **28**, 115-122.
40. Bernhart, S.H., Tafer, H., Mückstein, U., Flamm, C., Stadler, P.F., and Hofacker, I.L. (2006). Partition function and base pairing probabilities of RNA heterodimers. *Algorithms for Molecular Biology*, **1**, 1-10.
41. Cappannini, A., Ray, A., Purta, E., Mukherjee, S., Boccaletto, P., Moafinejad, S.N., Lechner, A., Barchet, C., Klaholz, B.P., and Stefaniak, F. (2024). MODOMICS: a database of RNA modifications and related information. 2023 update. *Nucleic acids research*, **52**, D239-D244.
42. Shen, T., Hu, Z., Sun, S., Liu, D., Wong, F., Wang, J., Chen, J., Wang, Y., Hong, L., and Xiao, J. (2024). Accurate RNA 3D structure prediction using a language model-based deep learning approach. *Nature methods*, 1-12.
43. Tsybulskyi, V., Semenchenko, E., and Meyer, I.M. (2023). e-RNA: a collection of web-servers for the prediction and visualisation of RNA secondary structure and their functional features. *Nucleic acids research*, **51**, W160-W167.
44. Gong, T., Ju, F., and Bu, D. (2024). Accurate prediction of RNA secondary structure including pseudoknots through solving minimum-cost flow with learned potentials. *Communications biology*, **7**, 297.
45. Hofacker, I.L., Fontana, W., Stadler, P.F., Bonhoeffer, L.S., Tacker, M., and Schuster, P. (1994). Fast folding and comparison of RNA secondary structures. *Monatshefte für Chemie*, **125**, 167-167.
46. Nawrocki, E.P., and Eddy, S.R. (2013). Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics*, **29**, 2933-2935.
47. Brázda, V., Kolomazník, J., Lýsek, J., Bartas, M., Fojta, M., Šťastný, J., and Mergny, J.-L. (2019). G4Hunter web application: a web server for G-quadruplex prediction. *Bioinformatics*, **35**, 3493-3495.
48. Rice, P., Longden, I., and Bleasby, A. (2000). EMBOSS: the European molecular biology open software suite. *Trends in genetics*, **16**, 276-277.
49. Macke, T.J., Ecker, D.J., Gutell, R.R., Gautheret, D., Case, D.A., and Sampath, R. (2001). RNAMotif, an RNA secondary structure definition and search algorithm. *Nucleic acids research*, **29**, 4724-4735.