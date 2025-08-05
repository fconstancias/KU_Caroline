# Introducing the Oral Cavity as a Novel Human In-Vivo Model System to Study Biofilm

## Repository Purpose

This repository is dedicated to centralizing the bioinformatics, data handling, and visualization code used in the project "Introducing the Oral Cavity as a Novel Human In-Vivo Model System to Study Biofilm." By bringing all related code into one place, we aim to streamline analysis workflows, ensure reproducibility, and facilitate collaboration among researchers involved in the project. This repository will contain scripts and conda environments for:

- Metagenomic data pre-processing:
	*  Quality trimming/ adaptor removal
	*  Filtering reads mapping the microbial masked version of the human genome

- Taxonomic and funcational profiling:
	*  MetaPhlAn/ HUMAnN profiling
	*  eHOMD taxonomic kraken/bracken profiling 
	*  	mOTUS taxonomic profiling 

- Assembly based appraoches
	*  co-assembly:
		*  megahit
		*  sample mapping with bowtie2: **check mappability**
		*  MMseq2 catalogue gene clustering
		*  Gene catalogue annotation: KEGG/COG/dbCAN

	*  single sample assembly: see <SqueezeMeta>
		*  megahit *vs* metaspades: test on few samples
		*  MMseq2 catalogue gene clustering: comparaison with co-assembly appraoch
		*  automatic binning: Metabat2
		*  dereplication with dRep
		*  MAGs taxonomy with GTDB
	
- Data visualization and Statistical analysis
- Result interpretation

We welcome contributions and collaborations to enhance and expand the codebase to further our understanding of biofilm dynamics in the oral cavity and beyond.

## Background

Biofilms play a central role in multiple oral and systemic diseases[^1]. A significant challenge in biofilm research is the accessibility of clinical samples, leading to a reliance on in vitro or animal models[^2]. To enhance prevention and treatment of biofilm-associated diseases, there is a need for human in-vivo models that can provide clinically translatable data.

Gingivitis, periodontitis, and dental caries are the most common biofilm-associated oral diseases[^3]. Inflammation drives gingivitis and periodontitis[^4], while frequent sugar intake contributes to dental caries[^5]. We have identified key characteristics of oral biofilm in periodontitis and dental caries compared to oral health[^6^][^7^][^8]. Additionally, we have studied biofilm compositional changes due to inadequate oral hygiene, frequent sugar intake, and periodontal treatment[^9^][^10^][^11]. However, the biological mechanisms of how oral biofilms adapt to changing environmental conditions remain unknown, including the bacterial genes regulating biofilm adaptation during disease-associated perturbations and post-treatment inflammation resolution.

Identifying bacterial genes that facilitate biofilm adaptation during inadequate oral hygiene, frequent sugar intake, and post-periodontal treatment is crucial for understanding biofilm responses to ecological changes in the oral cavity. These genes are potential molecular targets for preventing gingivitis, periodontitis, and dental caries. We hypothesize that baseline levels of specific bacterial gene signatures can be used for molecular risk assessment of oral diseases.

This project introduces the oral cavity as a novel in-vivo model system with various applications. It will demonstrate the feasibility of using bacterial genes as biomarkers for gingivitis, periodontitis, and dental caries, paving the way for molecular-based risk assessment and improved oral disease prevention. The project will also validate candidate biomarker genes in publicly available datasets to determine their relevance in non-oral inflammatory conditions, such as type 2 diabetes and cardiovascular diseases. This approach highlights the practicality of using the oral cavity as an in-vivo model system instead of in-vitro setups and animal models for biofilm-associated disease studies of both oral and non-oral origin.

## Aims
The overall aim is to demonstrate the feasibility of using the oral cavity as a human model system to study biofilms in-vivo. The specific aims are:

1. To identify bacterial genes affected during the initiation and resolution of oral inflammation and sugar stress.
2. To test if specific genes are biomarkers of different disease trajectories of gingivitis, periodontitis, and dental caries.
3. To validate the generalizability of identified genes in datasets of non-oral origin.

## Hypotheses
The general hypothesis is that the oral cavity is a unique but previously overlooked human model system for studying biofilms in-vivo. The specific hypotheses are:

1. Metagenomic datasets can identify bacterial genes involved in biofilm adaptation to oral inflammation and sugar stress.
2. Candidate genes associated with biofilm adaptation can be used as biomarkers to predict different disease trajectories of gingivitis, periodontitis, and dental caries.
3. Expression of candidate genes can be identified in biofilm retrieved from other body sites, such as the gut, as potential biomarkers for general medical disorders such as type 2 diabetes.

## Methods
### Sample Material
The project involves 360 saliva samples, 240 supragingival plaque samples, and 120 subgingival plaque samples, already collected as part of an ongoing PhD project. Samples were collected from 40 participants before (baseline), during (week 2), and after (week 4) 14 days of oral hygiene discontinuation and sugar stress. Additional samples were collected from 40 participants before (baseline), during (week 6), and after (week 12) non-surgical periodontal treatment. These samples provide an opportunity to identify bacterial genes regulating biofilm adaptation due to oral hygiene discontinuation, sugar stress, and non-surgical periodontal treatment.

### Metagenomic Sequencing
Data was generated at a collaborator's laboratory (ADM in Spain) using a protocol from previous studies[^12^][^13^][^14]. Shot-gun Illumina sequencing targeting a depth of 25 million read pairs was used to profile microbial genes down to low abundances[^15].

### Computational Biology and Statistical Approach
Functional profiling methods, including Metagenome Assembled Genomes (MAGs) binning, will be used to characterize the microbial gene composition[^12^][^13^][^14]. Data will identify genes up- or downregulated during inflammation initiation and resolution and sugar stress (Aim 1). Clinical endpoints (gingival inflammation at day 14 and week 6) will stratify responders in the three cohorts. Bacterial genes in baseline samples will test if gene signatures predict response trajectories to oral hygiene discontinuation, sugar stress, and non-surgical periodontal treatment (Aim 2). Publicly available datasets from stool, blood, and skin samples will validate bacterial genes associated with oral inflammation and sugar stress in general medical conditions such as type 2 diabetes and cardiovascular disease (Aim 3).

## Practical Circumstances and the Role of the Candidate
The PhD candidate will be responsible for data analysis in accordance with Aims 1-3, working closely with two current PhD students: Vincent Vendius (bioinformatician) and Christine Lundtorp (sample collection and clinical interpretation). The project will be conducted within an interdisciplinary research group with expertise in bioinformatics, molecular biology, and odontology, led by Professor Daniel Belstrøm (main supervisor) and co-supervisors Professors Michael Givskov, Martin Sikora, and Hannes Schrøder. External collaborators include Professor Sebastian Schlafer and bioinformatician Florentin Constancias.

All samples have been collected, and molecular analyses completed. Salary costs are covered by an internal grant from the Department of Odontology, University of Copenhagen. Sample collection and metagenomic analysis costs are funded by Innovation Fund Denmark and ADM, ensuring no additional operational costs.

## Ethical Approval and Data Management
The project has been approved by the regional ethical committee of the capital region of Denmark (H-21003295) and reported to the local data authority of SUND (514-0649/21-3000). It is registered in clinical trials.gov (NCT05073393, NCT05268757, NCT05518747). Collaboration and data handling agreements have been signed with all parties involved.

## Publication Strategy
Data from this PhD project is expected to be published in at least three international peer-reviewed journal papers. Papers based on Aims 1-2 will be submitted to high-impact journals in oral sciences, such as the Journal of Dental Research and the Clinical Journal of Periodontology. The paper addressing the generalizability of gene expression (Aim 3) will target a broad microbiology journal, such as Nature Microbiology. The PhD candidate will be the first author on all papers. Additionally, the PhD student will present data at national and international conferences, such as the annual meetings of the International Association of Dental Research (IADR). Conference attendance and open access fees are covered by the Department of Odontology, University of Copenhagen.

## Risk Assessment and Alternative Approaches
All data has been generated and is available to the PhD student, who will work within a supportive environment of experienced supervisors and fellow PhD students. Hence, the current project is considered low risk, with no need for alternative approaches.

## Perspectives
This project uses the oral cavity as a human in-vivo model system to identify bacterial genes involved in biofilm adaptation during early and chronic oral inflammation and sugar stress. By focusing on genes rather than bacterial species, it will be possible to compare oral data with non-oral data. Validation in non-oral datasets will demonstrate that findings from the oral cavity can translate into non-oral conditions. This project has the potential to revolutionize risk assessment in dental clinics, laying the foundation for saliva-based molecular risk assessment of oral diseases at preclinical stages. Such advancements will improve prevention efforts and have a substantial societal impact by better preventing oral diseases.

## Time Schedule

<!-- Insert your time schedule here -->

## References
[^1]: Vestby LK, Grønseth T, Simm R, Nesse LL. Bacterial Biofilm and its Role in the Pathogenesis of Disease. Antibiotics (Basel). 2020 Feb 3;9(2):59.
[^2]: Jakobsen TH, Bjarnsholt T, Jensen PØ, Givskov M, Høiby N. Targeting quorum sensing in Pseudomonas aeruginosa biofilms: current and emerging inhibitors. Future Microbiol. 2013 Jul;8(7):901-21. 
[^3]: Marsh PD, Zaura E. Dental biofilm: ecological interactions in health and disease. J Clin Periodontol. 2017 Mar;44 Suppl 18:S12-S22.
[^4]: Van Dyke TE, Bartold PM, Reynolds EC. The Nexus Between Periodontal Inflammation and Dysbiosis. Front Immunol. 2020 Mar 31;11:511.
[^5]: Selwitz RH, Ismail AI, Pitts NB. Dental caries. Lancet. 2007 Jan 6;369(9555):51-9. 
[^6]: Belstrøm D, Fiehn NE, Nielsen CH, Kirkby N, Twetman S, Klepac-Ceraj V, Paster BJ, Holmstrup P. Differences in bacterial saliva profile between periodontitis patients and a control cohort. J Clin Periodontol. 2014 Feb;41(2):104-12.
[^7]: Belstrøm D, Fiehn NE, Nielsen CH, Holmstrup P, Kirkby N, Klepac-Ceraj V, Paster BJ, Twetman S. Altered bacterial profiles in saliva from adults with caries lesions: a case-cohort study. Caries Res. 2014;48(5):368-75.
[^8]: Belstrøm D, Paster BJ, Fiehn NE, Bardow A, Holmstrup P. Salivary bacterial fingerprints of established oral disease revealed by the Human Oral Microbe Identification using Next Generation Sequencing (HOMINGS) technique. J Oral Microbiol. 2016 Jan 14;8:30170.
[^9]: Belstrøm D, Grande MA, Sembler-Møller ML, Kirkby N, Cotton SL, Paster BJ, Holmstrup P. Influence of periodontal treatment on subgingival and salivary microbiotas. J Periodontol. 2018 May;89(5):531-539.
[^10]: Belstrøm D, Sembler-Møller ML, Grande MA, Kirkby N, Cotton SL, Paster BJ, Twetman S, Holmstrup P. Impact of Oral Hygiene Discontinuation on Supragingival and Salivary Microbiomes. JDR Clin Trans Res. 2018 Jan;3(1):57-64.
[^11]: Lundtorp-Olsen C, Enevold C, Juel Jensen CA, Stofberg SN, Twetman S, Belstrøm D. Impact of Probiotics on the Salivary Microbiota and Salivary Levels of Inflammation-Related Proteins during Short-Term Sugar Stress: A Randomized Controlled Trial. Pathogens. 2021 Mar 25;10(4):392. 
[^12]: Belstrøm D, Constancias F, Liu Y, Yang L, Drautz-Moses DI, Schuster SC, Kohli GS, Jakobsen TH, Holmstrup P, Givskov M. Metagenomic and metatranscriptomic analysis of saliva reveals disease-associated microbiota in patients with periodontitis and dental caries. NPJ Biofilms Microbiomes. 2017 Oct 2;3:23.
[^13]: Belstrøm D, Constancias F, Drautz-Moses DI, Schuster SC, Veleba M, Mahé F, Givskov M. Periodontitis associates with species-specific gene expression of the oral microbiota. NPJ Biofilms Microbiomes. 2021 Sep 23;7(1):76.
[^14]: Belstrøm D, Constancias F, Markvart M, Sikora M, Sørensen CE, Givskov M. Transcriptional Activity of Predominant Streptococcus Species at Multiple Oral Sites Associate With Periodontal Status. Front Cell Infect Microbiol. 2021 Sep 21;11:752664.
[^15]: Pereira-Marques J, Hout A, Ferreira RM, Weber M, Pinto-Ribeiro I, van Doorn LJ, et al. Impact of Host DNA and Sequencing Depth on the Taxonomic Resolution of Whole Metagenome Sequencing for Microbiome Analysis. Front Microbiol. 2019;10:1277.
[^16]: Turnbaugh PJ, Ley RE, Hamady M, Fraser-Liggett CM, Knight R, Gordon JI. The human microbiome project. Nature. 2007 Oct 18;449(7164):804-10.
[^17]: Perez-Riverol Y, Csordas A, Bai J, Bernal-Llinares M, Hewapathirana S, Kundu DJ, Inuganti A, Griss J, Mayer G, Eisenacher M, Pérez E, Uszkoreit J, Pfeuffer J, Sachsenberg T, Yilmaz S, Tiwary S, Cox J, Audain E, Walzer M, Jarnuczak AF, Ternent T, Brazma A, Vizcaíno JA. The PRIDE database and related tools and resources in 2019: improving support for quantification data. Nucleic Acids Res. 2019 Jan 8;47(D1):D442-D450.
