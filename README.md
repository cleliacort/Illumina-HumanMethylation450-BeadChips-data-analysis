# Illumina-HumanMethylation450-BeadChips-data-analysis

The aim of this project is to provide a guide for analyze data coming from the Illumina Infinium HumanMethylation 450k BeadChip, using specific R package. 

The 450k BeadChip contains about 485512 probes covering gene region and CpG island, plus additional informative site (CpG sites outside of CpG island, Non-CpG methylated sites identifies in human stem cells). To have a more compressive view of the methylation state, the gene coverage was targeted for sites across the different gene regions like the promoter region, the 5’ UTR, gene body and 3’ UTR.

The physical support is a silicon slide and in our case it contains 8 arrays named on the basis of their specific position. Each array corresponds to one sample. Before doing the hybridization, the DNA samples are treated with the bisulphite to maintain the methylated cytosines once it is amplified. After this treatment,the bisulphite-converted DNA is amplify and then it is fragmented with some enzyme. 
Once the previous steps were executed, the bisulphite-converted DNA fragments are put on the slide for making the hybridization with the probes. A specific scan were used to detect the hybridization between the sample and the probes. It was able to measure the signal
intensity released from each well. The 450k BeadChip is based on two types of chemistry which are Infinium I and Infinium II. They differ in the types of bead, the number of probes and the number of colors used to detect the hybridization.
