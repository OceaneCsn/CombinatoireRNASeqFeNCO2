library(CoRegNet)
library(stringr)

data(CIT_BLCA_EXP,HumanTF,CIT_BLCA_Subgroup)
# expression matrix
dim(CIT_BLCA_EXP)
head(intersect(rownames(CIT_BLCA_EXP),HumanTF))

# construction du réseau de coregulation
# mettre ici notre discretisation perso si besoin
options("mc.cores"=6)
grn = hLICORN(head(CIT_BLCA_EXP,100), TFlist=HumanTF, parallel = "multicore", verbose = T)
print(grn)

# Uses a network in the form of a coregnet object to compute regulatory influence to estimate 
#the transcriptional activity of each regulators in each sample of the given expression data.
# addEvidences to add external evidences
influence = regulatorInfluence(grn,head(CIT_BLCA_EXP,200))


# Based on the frequency and specificity of co-regulation, this functions extracts from a coregnet 
# network all the cooperative regulators.
coregs= coregulators(grn)


display(grn,expressionData=head(CIT_BLCA_EXP,200),TFA=influence)


# On essaie avec nos données hihi, application aux gènes qui répondent au CO2

load("normalized.count_At.RData")
load("GenesNitrate_At.RData")

TF <- read.table("TFs_PlnTFDB.txt", h=T, sep = '\t')
TF$AGI <- str_split_fixed(TF$Protein.ID, '\\.', 2)[,1]
data <- normalized.count[sharedBy3,]
intersect(rownames(normalized.count[sharedBy3,]),as.vector(TF$AGI))
options("mc.cores"=6)

discexp = discretizeExpressionData(normalized.count[sharedBy3,],refSamples=c("cNF_1", "cNF_2", "cNF_3"))


grn = hLICORN(normalized.count[sharedBy3,], TFlist=TF$AGI, parallel = "multicore", verbose = T, minGeneSupport = 0)
print(grn)
