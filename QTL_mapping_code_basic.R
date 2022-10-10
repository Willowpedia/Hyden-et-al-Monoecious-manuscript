library(qtl)
#Area
area<-read.cross("csv", file='S2_8_area.csv', genotypes=c("AA", "AB", "BB"), alleles=c("A", "B"), estimate.map=TRUE)
area_MQM<-sim.geno(area, step=2, n.draws=100, error.prob=0.01)
area_MQM_1<-scanone(area_MQM, method="imp") #QTL plot, no covariate
operm_area_1<-scanone(area_MQM, method="imp", n.perm=1000) #find number of permutations for threshold, this takes a while
summary(operm_area_1)
#choose value that is 5% threshold
summary(operm_area_1)[1]
#3.89
#plot without covariate
pdf('area.pdf', width=24, height=6)
plot(area_MQM_1)
abline(summary(operm_area_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()
write.csv(area_MQM_1, "area_results.csv")

#Height
height<-read.cross("csv", file='S2_8_height.csv', genotypes=c("AA", "AB", "BB"), alleles=c("A", "B"), estimate.map=TRUE)
height_MQM<-sim.geno(height, step=2, n.draws=100, error.prob=0.01)
height_MQM_1<-scanone(height_MQM, method="imp") #QTL plot, no covariate
operm_height_1<-scanone(height_MQM, method="imp", n.perm=1000) #find number of permutations for threshold, this takes a while
summary(operm_height_1)
#choose value that is 5% threshold
summary(operm_height_1)[1]
#3.86
#plot without covariate
pdf('height.pdf', width=24, height=6)
plot(height_MQM_1)
abline(summary(operm_height_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()
write.csv(height_MQM_1, "height_results.csv")

#SPAD
spad<-read.cross("csv", file='S2_8_spad.csv', genotypes=c("AA", "AB", "BB"), alleles=c("A", "B"), estimate.map=TRUE)
spad_MQM<-sim.geno(spad, step=2, n.draws=100, error.prob=0.01)
spad_MQM_1<-scanone(spad_MQM, method="imp") #QTL plot, no covariate
operm_spad_1<-scanone(spad_MQM, method="imp", n.perm=1000) #find number of permutations for threshold, this takes a while
summary(operm_spad_1)
#choose value that is 5% threshold
summary(operm_spad_1)[1]
#7.98
#plot without covariate
pdf('spad.pdf', width=24, height=6)
plot(spad_MQM_1)
abline(summary(operm_spad_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()
write.csv(spad_MQM_1, "spad_results.csv")

#Stem Color
stemcolor<-read.cross("csv", file='S2_8_stemcolor.csv', genotypes=c("AA", "AB", "BB"), alleles=c("A", "B"), estimate.map=TRUE)
stemcolor_MQM<-sim.geno(stemcolor, step=2, n.draws=100, error.prob=0.01)
stemcolor_MQM_1<-scanone(stemcolor_MQM, method="imp") #QTL plot, no covariate
operm_stemcolor_1<-scanone(stemcolor_MQM, method="imp", n.perm=1000) #find number of permutations for threshold, this takes a while
summary(operm_stemcolor_1)
#choose value that is 5% threshold
summary(operm_stemcolor_1)[1]
#4.18
#plot without covariate
pdf('stemcolor.pdf', width=24, height=6)
plot(stemcolor_MQM_1)
abline(summary(operm_stemcolor_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()
write.csv(stemcolor_MQM_1, "stemcolor_results.csv")

#Survival
survival<-read.cross("csv", file='S2_8_survival.csv', genotypes=c("AA", "AB", "BB"), alleles=c("A", "B"), estimate.map=TRUE)
survival_MQM<-sim.geno(survival, step=2, n.draws=100, error.prob=0.01)
survival_MQM_1<-scanone(survival_MQM, method="imp") #QTL plot, no covariate
operm_survival_1<-scanone(survival_MQM, method="imp", n.perm=1000) #find number of permutations for threshold, this takes a while
summary(operm_survival_1)
#choose value that is 5% threshold
summary(operm_survival_1)[1]
#6.96
#plot without covariate
pdf('survival.pdf', width=24, height=6)
plot(survival_MQM_1)
abline(summary(operm_survival_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()
write.csv(survival_MQM_1, "survival_results.csv")

#Sex
sex<-read.cross("csv", file='S2_8_sex.csv', genotypes=c("AA", "AB", "BB"), alleles=c("A", "B"), estimate.map=TRUE)
sex_MQM<-sim.geno(sex, step=2, n.draws=100, error.prob=0.01)
sex_MQM_1<-scanone(sex_MQM, method="imp", model="binary") #QTL plot, no covariate
operm_sex_1<-scanone(sex_MQM, method="imp", n.perm=1000) #find number of permutations for threshold, this takes a while
summary(operm_sex_1)
#choose value that is 5% threshold
summary(operm_sex_1)[1]
#3.83
#plot without covariate
pdf('sex.pdf', width=24, height=6)
plot(sex_MQM_1)
abline(summary(operm_sex_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()
write.csv(sex_MQM_1, "sex_results.csv")

max(sex_MQM_1) #find maxiumum marker for covariate
g<-pull.geno(fill.geno(sex))[, "S15_6318849"] #pull covariate
sex_MQM.c15<-scanone(sex_MQM, model="binary", addcovar=g) #QTL with covaraite
operm_sex<-scanone(sex_MQM, model="binary", addcovar=g, n.perm=1000) #find number of permutations for threshold, this takes a while
plot(operm_sex) #visualize distributions of permutation tests
summary(operm_sex)
#choose value that is 5% threshold
summary(operm_sex)[1]
#3.62
#plot with covariate
pdf('sex_MQM.c15.pdf', width=24, height=6)
plot(sex_MQM.c15)
abline(summary(operm_sex)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()
write.csv(sex_MQM.c15, "sex_covariate_results.csv")



par(mfrow=c(5,1))
pdf('QTL_plots_no_covariate_2022.pdf', width=24, height=6)
plot(area_MQM_1)
abline(summary(operm_area_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
plot(height_MQM_1)
abline(summary(operm_height_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
plot(spad_MQM_1)
abline(summary(operm_spad_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
plot(stemcolor_MQM_1)
abline(summary(operm_stemcolor_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
plot(survival_MQM_1)
abline(summary(operm_survival_1)[1],0, col='red') #abline is summary of operm, this is the permutation-based threshold
dev.off()

