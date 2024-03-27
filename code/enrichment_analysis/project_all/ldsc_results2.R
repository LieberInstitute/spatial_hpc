

# read data 
dat <- read.table("ldsc_results.txt",as.is=T,header=T,sep="\t")

traits <- c(
"ADHD",
"Alzheimer Disease3",
 "Anorexia",
 "Autism",
 "BMI",
 "Bipolar Disorder2",
 "GSCAN_AgeSmk",
 "GSCAN_CigDay",
 "GSCAN_DrnkWk",
 "GSCAN_SmkCes",
 "GSCAN_SmkInit",
 "Height",
 "Type_2_Diabetes",
 "mdd2019edinburgh",
 "epilepsy",
 "Intelligence",
 "Parkinson Disease",
 "Schizophrenia_PGC3",
 "Education Years",
 "Neuroticism"
)
idx <- is.element(dat$trait,traits)
dat <- dat[idx,]

# rename traits
dat$trait[dat$trait=="GSCAN_AgeSmk"] <- "Age of smoking"
dat$trait[dat$trait=="GSCAN_CigDay"] <- "Cigarettes per day"
dat$trait[dat$trait=="GSCAN_DrnkWk"] <- "Drinks per week"
dat$trait[dat$trait=="GSCAN_SmkCes"] <- "Smoking cessation"
dat$trait[dat$trait=="GSCAN_SmkInit"] <- "Smoking initiation"
dat$trait[dat$trait=="Schizophrenia_PGC3"] <- "Schizophrenia"
dat$trait[dat$trait=="Bipolar Disorder2"] <- "Bipolar"
dat$trait[dat$trait=="epilepsy"] <- "Epilepsy"
dat$trait[dat$trait=="Type_2_Diabetes"] <- "Type 2 Diabetes"
dat$trait[dat$trait=="Alzheimer Disease3"] <- "Alzheimer Disease"
dat$trait[dat$trait=="mdd2019edinburgh"] <- "Depression"

# FDR
dat$p_zcore <- pnorm(abs(dat$Coefficient_z.score),lower.tail=F)*2
dat$FDR <- p.adjust(dat$p_zcore,method="fdr")


write.csv(dat,"ldsc_results.csv")






