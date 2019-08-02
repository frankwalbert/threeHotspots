library(tidyverse)
library(magrittr)
library(tidyselect)
library(gridExtra)
library(readxl)

# this file has the non-detected fatty acids removed
# can create this from Supplementary data in paper
dat <- read_excel("fattyAcids_reformat_190328_deleteNACols.xlsx", col_names = TRUE)

# start with NEFAs
# make the NA columns numeric
dat %<>% mutate_at(vars(ends_with("_nmol")), funs(as.numeric(.)))
# make a total NEFA column
dat %<>% mutate(
  totalNEFA_nmol = rowSums(select(., ends_with("_nmol")), na.rm=TRUE),
  totalSaturatedNEFA_nmol = myristic_14_0_nmol + palmitic_16_0_nmol + stearic_18_0_nmol,
  totalUnsaturatedNEFA_nmol = palmitoleic_16_1_nmol + oleic_18_1_nmol,
  totalSaturated_percent = myristic_14_0_percent + palmitic_16_0_percent + stearic_18_0_percent,
  totalUnsaturated_percent = palmitoleic_16_1_percent + oleic_18_1_percent,
  total16_nmol = palmitic_16_0_nmol + palmitoleic_16_1_nmol,
  total18_nmol = stearic_18_0_nmol + oleic_18_1_nmol,
  total16_percent = palmitic_16_0_percent + palmitoleic_16_1_percent,
  total18_percent = stearic_18_0_percent + oleic_18_1_percent
)
dat <- as.data.frame(dat)
dat$grouping <- as.factor(dat$grouping)
# change factor order for the genotypes
dat$genotype <- factor(dat$genotype, levels=c("BY", "BY_OAF1_L63S", "BY_OLE1_FAR_RM", "BY_OAF1_L63S_OLE1_FAR_RM" , "RM"))

testGenotypes <- c("RM", "BY_OAF1_L63S","BY_OLE1_FAR_RM","BY_OAF1_L63S_OLE1_FAR_RM")
testMeasurements <- c(vars_select(names(dat), ends_with("_nmol")), "protein_ug", "OD")
names(testMeasurements)[(length(testMeasurements) - 1):length(testMeasurements)] <- c("protein_ug", "OD")

testMeasurementsComposition <- c(vars_select(names(dat), ends_with("_percent")))


# residuals with means added back in so that plots are more interpretable
datResidualNEFAs <- dat %>% mutate(
    myristic_14_0_nmol = mean(myristic_14_0_nmol) + lm(myristic_14_0_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    palmitoleic_16_1_nmol = mean(palmitoleic_16_1_nmol) + lm(palmitoleic_16_1_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    palmitic_16_0_nmol = mean(palmitic_16_0_nmol) + lm(palmitic_16_0_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    oleic_18_1_nmol = mean(oleic_18_1_nmol) + lm(oleic_18_1_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    stearic_18_0_nmol = mean(stearic_18_0_nmol) + lm(stearic_18_0_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    totalSaturatedNEFA_nmol = mean(totalSaturatedNEFA_nmol) + lm(totalSaturatedNEFA_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    totalUnsaturatedNEFA_nmol = mean(totalUnsaturatedNEFA_nmol) + lm(totalUnsaturatedNEFA_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    total16_nmol = mean(total16_nmol) + lm(total16_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    total18_nmol = mean(total18_nmol) + lm(total18_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    totalNEFA_nmol = mean(totalNEFA_nmol) + lm(totalNEFA_nmol ~ grouping + acquisitionOrder_NEFA + protein_ug)$residuals,
    protein_ug = mean(protein_ug) + lm(protein_ug ~ grouping + acquisitionOrder_NEFA)$residuals
)


# make a plot for each component, showing the different genotypes
pdf("individualTraitsResiduals.pdf", width=15, height=15)
#ggPlots <- lapply(testMeasurements[c(7, 8, 9, 10, 3, 2, 5, 4, 1, 6, 11, 12)], function(i){
ggPlots <- lapply(testMeasurements, function(i){
  #ggplot(dat, aes(x=genotype, y=eval(parse(text=i)))) +
  ggplot(datResidualNEFAs, aes(x=genotype, y=eval(parse(text=i)))) +
    geom_jitter(position = position_jitter(width=0.05)) +
    stat_summary(fun.y = mean, fun.ymin=mean, fun.ymax=mean, colour = "red", geom="crossbar", width=0.5) +
    ylab(i) + xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
})
grid.arrange(ggPlots[[1]], ggPlots[[2]], ggPlots[[3]], ggPlots[[4]], ggPlots[[5]], ggPlots[[6]], ggPlots[[7]], ggPlots[[8]], ggPlots[[9]], ggPlots[[10]], nrow=3, ncol=4)
dev.off()


sapply(testGenotypes, function(otherGenotype){
    sapply(testMeasurements[c(7, 8, 3, 2, 5, 4, 1, 6, 9, 10)], function(thisMeasurement){
      ret <- NA
      
      # run with the residuals:
      group1 <- datResidualNEFAs[datResidualNEFAs$genotype=="BY", thisMeasurement]
      group2 <- datResidualNEFAs[datResidualNEFAs$genotype==otherGenotype, thisMeasurement]

      if(length(which(!is.na(group1))) >= 3 & length(which(!is.na(group2))) >= 3){
          #pVals
        ret <- t.test(group1, group2)$p.value
        # or estimates
        ret <- t.test(group1, group2)$est[2] - t.test(group1, group2)$est[1]
    }
      ret
  })
})

# report the RM column in the paper; other genotypes are better modeled with an ANOVA (see below)
#                         RM            BY_OAF1_L63S    BY_OLE1_FAR_RM        BY_OAF1_L63S_OLE1_FAR_RM
#totalSaturatedNEFA_nmol   1.607057e-02   0.58789491   0.0011104113              0.002204363
#totalUnsaturatedNEFA_nmol 1.210468e-06   0.58532268   0.1180399913              0.211196408
#palmitic_16_0_nmol        1.609060e-01   0.48363332   0.0020818167              0.004129041
#palmitoleic_16_1_nmol     2.625278e-05   0.97757199   0.0945852911              0.427231953
#stearic_18_0_nmol         9.300606e-01   0.15408275   0.0002536384              0.001172428
#oleic_18_1_nmol           3.333311e-09   0.02386221   0.5637599673              0.027022409
#myristic_14_0_nmol        1.020857e-04   0.66740382   0.0009036808              0.023591579
#totalNEFA_nmol            9.618905e-07   0.72863837   0.0082467757              0.028510678
#total16_nmol              1.506907e-05   0.83421238   0.0052638599              0.044585566
#total18_nmol              1.126304e-08   0.01839508   0.2472087116              0.016329592

# effect size estimates
#                                      RM           BY_OAF1_L63S    BY_OLE1_FAR_RM          BY_OAF1_L63S_OLE1_FAR_RM
#totalSaturatedNEFA_nmol.mean of y    0.914865670  -0.16400828      2.2415439                1.7491470
#totalUnsaturatedNEFA_nmol.mean of y 13.147927138   0.55260954      1.6049497                1.2755483
#palmitic_16_0_nmol.mean of y         0.420123841  -0.21707020      1.9389587                1.5060429
#palmitoleic_16_1_nmol.mean of y      7.542083544   0.02431227      1.4789884                0.6352971
#stearic_18_0_nmol.mean of y          0.003770833   0.06888622      0.1395066                0.1250814
#oleic_18_1_nmol.mean of y            5.605843594   0.52829726      0.1259613                0.6402512
#myristic_14_0_nmol.mean of y         0.490970996  -0.01582431      0.1630786                0.1180227
#totalNEFA_nmol.mean of y            14.062792808   0.38860126      3.8464935                3.0246953
#total16_nmol.mean of y               7.962207384  -0.19275792      3.4179472                2.1413400
#total18_nmol.mean of y               5.609614428   0.59718349      0.2654678                0.7653326


# model the double-swap as an interaction:
#lapply(testMeasurements[c(7, 8, 3, 2, 5, 4, 1, 6, 9, 10)], function(thisMeasurement){
NEFA_pVals <- t(sapply(testMeasurements, function(thisMeasurement){
    thisDat <- datResidualNEFAs[datResidualNEFAs$genotype != "RM",]
    #thisDat <- dat[dat$genotype != "RM", c("OAF1", "OLE1", vars_select(names(dat), ends_with("_nmol")))]
    pheno <- thisDat[, thisMeasurement]
    ret <- lm(pheno ~ thisDat$OAF1 * thisDat$OLE1)
    anova(ret)[1:3,5]
}))
colnames(NEFA_pVals) <- c("OAF1", "OLE1", "interaction")
#                           OAF1         OLE1         interaction
#myristic_14_0_nmol        0.263337695 3.589160e-05   0.5854731
#palmitoleic_16_1_nmol     0.436829112 5.885573e-02   0.4106554
#palmitic_16_0_nmol        0.238518236 3.582241e-06   0.6896880
#oleic_18_1_nmol           0.004788591 4.660302e-01   0.9654749
#stearic_18_0_nmol         0.271257275 8.457361e-04   0.1004838
#totalNEFA_nmol            0.787686454 8.394062e-04   0.4552745
#totalSaturatedNEFA_nmol   0.247434748 1.074948e-06   0.5565400
#totalUnsaturatedNEFA_nmol 0.861650353 8.331829e-02   0.4940875
#total16_nmol              0.263315536 3.357261e-04   0.4050719
#total18_nmol              0.004764330 2.138839e-01   0.7751521
#protein_ug                0.057573245 6.705865e-01   0.3777256
#OD                        0.041172166 4.811750e-03   0.3875861


NEFA_coefs <- t(sapply(testMeasurements, function(thisMeasurement){
    thisDat <- datResidualNEFAs[datResidualNEFAs$genotype != "RM",]
    #thisDat <- dat[dat$genotype != "RM", c("OAF1", "OLE1", vars_select(names(dat), ends_with("_nmol")))]
    pheno <- thisDat[, thisMeasurement]
    ret <- lm(pheno ~ thisDat$OAF1 * thisDat$OLE1)
    ret$coef
}))
# if we want to see coefficients (for effect directions)
#                           (Intercept) thisDat$OAF1RM thisDat$OLE1RM thisDat$OAF1RM:thisDat$OLE1RM
#myristic_14_0_nmol          0.4675433    -0.01582431      0.1630786                   -0.02923158
#palmitoleic_16_1_nmol      13.0686827     0.02431227      1.4789884                   -0.86800361
#palmitic_16_0_nmol          5.9612536    -0.21707020      1.9389587                   -0.21584562
#oleic_18_1_nmol             3.8356347     0.52829726      0.1259613                   -0.01400732
#stearic_18_0_nmol           0.6739285     0.06888622      0.1395066                   -0.08331135
#totalNEFA_nmol             24.0070428     0.38860126      3.8464935                   -1.21039949
#totalSaturatedNEFA_nmol     7.1027254    -0.16400828      2.2415439                   -0.32838856
#totalUnsaturatedNEFA_nmol  16.9043174     0.55260954      1.6049497                   -0.88201093
#total16_nmol               19.0299363    -0.19275792      3.4179472                   -1.08384923
#total18_nmol                4.5095632     0.59718349      0.2654678                   -0.09731868
#protein_ug                394.6839743  -100.74830743    -45.7344795                   61.90462587
#OD                          0.3270000     0.01540000      0.0206000                   -0.00880000



# plots
# have to reshape all the nmol values
datGather <- gather(datResidualNEFAs[,c("sample", "genotype", vars_select(names(datResidualNEFAs), ends_with("_nmol")))], names(vars_select(names(datResidualNEFAs), ends_with("_nmol"))), key="compound", value="nmol")
# change factor order
datGather$compound <- factor(datGather$compound, levels=c("myristic_14_0_nmol", "palmitic_16_0_nmol", "palmitoleic_16_1_nmol", "stearic_18_0_nmol", "oleic_18_1_nmol", "total16_nmol", "total18_nmol", "totalSaturatedNEFA_nmol", "totalUnsaturatedNEFA_nmol", "totalNEFA_nmol"))

# plot!
pdf("allCombined_NEFA.pdf", width=11, height=6)
print(ggplot(datGather[datGather$genotype !="RM" & !str_detect(datGather$compound, "total"),], aes(x=compound, y=nmol, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

print(ggplot(datGather[!str_detect(datGather$compound, "total"),], aes(x=compound, y=nmol, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF", "purple")) + theme_light())

print(ggplot(datGather[datGather$genotype !="RM" & str_detect(datGather$compound, "total"),], aes(x=compound, y=nmol, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

print(ggplot(datGather[str_detect(datGather$compound, "total"),], aes(x=compound, y=nmol, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF", "purple")) + theme_light())

print(ggplot(datGather[datGather$genotype !="RM" & str_detect(datGather$compound, "total16|total18"),], aes(x=compound, y=nmol, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

print(ggplot(datGather[datGather$genotype !="RM" & str_detect(datGather$compound, "totalSaturated|totalUnsaturated"),], aes(x=compound, y=nmol, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

print(ggplot(datGather[datGather$genotype !="RM" & str_detect(datGather$compound, "totalNEFA_nmol"),], aes(x=compound, y=nmol, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

dev.off()







####################################
# lipid composition

# raw
datGather <- gather(dat[,c("sample", "genotype", vars_select(names(dat), ends_with("_percent")))], names(vars_select(names(dat), ends_with("_percent"))), key="compound", value="fraction")
# change factor order
datGather$compound <- factor(datGather$compound, levels=c("myristic_14_0_percent", "palmitic_16_0_percent", "palmitoleic_16_1_percent", "stearic_18_0_percent", "oleic_18_1_percent", "total16_percent", "total18_percent", "totalSaturated_percent", "totalUnsaturated_percent"))
# remove samples with no oleic acid (which isn't believable)
datGather <- filter(datGather, !sample %in% c("3_3", "5_5", "4_5"))

# residuals with means added for easier interpration:
datResidualComposition <- dat %>% mutate(
    myristic_14_0_percent = mean(myristic_14_0_percent) + lm(myristic_14_0_percent ~ grouping + acquisitionOrder_composition)$residuals,
    palmitoleic_16_1_percent = mean(palmitoleic_16_1_percent) + lm(palmitoleic_16_1_percent ~ grouping + acquisitionOrder_composition)$residuals,
    palmitic_16_0_percent = mean(palmitic_16_0_percent) + lm(palmitic_16_0_percent ~ grouping + acquisitionOrder_composition)$residuals,
    oleic_18_1_percent = mean(oleic_18_1_percent) + lm(oleic_18_1_percent ~ grouping + acquisitionOrder_composition)$residuals,
    stearic_18_0_percent = mean(stearic_18_0_percent) + lm(stearic_18_0_percent ~ grouping + acquisitionOrder_composition)$residuals,
    total16_percent = mean(total16_percent) + lm(total16_percent ~ grouping + acquisitionOrder_composition)$residuals,
    total18_percent = mean(total18_percent) + lm(total18_percent ~ grouping + acquisitionOrder_composition)$residuals,
    totalSaturated_percent = mean(totalSaturated_percent) + lm(totalSaturated_percent ~ grouping + acquisitionOrder_composition)$residuals,
    totalUnsaturated_percent = mean(totalUnsaturated_percent) + lm(totalUnsaturated_percent ~ grouping + acquisitionOrder_composition)$residuals
)

# statistics
# let's just treat the summed un/saturated fatty acid fraction as a numeric value and do ANOVA/T-test
composition_pVals <- t(sapply(testMeasurementsComposition, function(thisMeasurement){
  thisDat <- datResidualComposition[datResidualComposition$genotype != "RM",]
  #thisDat <- dat[dat$genotype != "RM",]
  pheno <- thisDat[, thisMeasurement]
  ret <- lm(pheno ~ thisDat$OAF1 * thisDat$OLE1)
  anova(ret)[1:3,5]
}))
colnames(composition_pVals) <- c("OAF1", "OLE1", "interaction")
#                         OAF1        OLE1        interaction
#myristic_14_0_percent    0.391923764 0.9933897   0.1544754
#palmitoleic_16_1_percent 0.182051107 0.4895209   0.5062466
#palmitic_16_0_percent    0.798767378 0.4453890   0.2882074
#oleic_18_1_percent       0.313152013 0.9908536   0.1506852
#stearic_18_0_percent     0.327848065 0.8894291   0.2389349
#totalSaturated_percent   0.811814501 0.7119834   0.2487220
#totalUnsaturated_percent 0.773116315 0.6866427   0.2603495
#total16_percent          0.001598449 0.9328807   0.4113122
#total18_percent          0.001339001 0.8293150   0.3339780

# coefficients
composition_coefs <- t(sapply(testMeasurementsComposition, function(thisMeasurement){
    thisDat <- datResidualComposition[datResidualComposition$genotype != "RM",]
    #thisDat <- dat[dat$genotype != "RM",]
    pheno <- thisDat[, thisMeasurement]
    ret <- lm(pheno ~ thisDat$OAF1 * thisDat$OLE1)
    ret$coef
}))
colnames(composition_coefs) <- c("intercept", "OAF1", "OLE1", "interaction")
#intercept         OAF1         OLE1 interaction
#myristic_14_0_percent    0.02104234  0.001759857  0.004255007 -0.00855820
#palmitoleic_16_1_percent 0.32017921 -0.076406825 -0.051083133  0.05007456
#palmitic_16_0_percent    0.29218634  0.029403944  0.065894629 -0.07696997
#oleic_18_1_percent       0.19736072 -0.014848193 -0.047522073  0.09578304
#stearic_18_0_percent     0.16942957  0.059327079  0.028752835 -0.06501305
#totalSaturated_percent   0.48265825  0.090490880  0.098902471 -0.15054121
#totalUnsaturated_percent 0.51753993 -0.091255018 -0.098605206  0.14585759
#total16_percent          0.61236556 -0.047002881  0.014811496 -0.02689541
#total18_percent          0.36679028  0.044478886 -0.018769238  0.03076999

# to get the BY vs RM comparison
sapply(testGenotypes, function(otherGenotype){
  sapply(testMeasurementsComposition, function(thisMeasurement){
    ret <- NA
    group1 <- datResidualComposition[datResidualComposition$genotype=="BY", thisMeasurement]
    group2 <- datResidualComposition[datResidualComposition$genotype==otherGenotype, thisMeasurement]
    
    if(length(which(!is.na(group1))) >= 3 & length(which(!is.na(group2))) >= 3){
        ret <- t.test(group1, group2)$p.value
        # or get the estimate (other genotype minus BY)
        ret <- t.test(group1, group2)$est[2] - t.test(group1, group2)$est[1]
        print(t.test(group1, group2))
    }
    ret
  })
})
# pVals
#                             RM            BY_OAF1_L63S   BY_OLE1_FAR_RM         BY_OAF1_L63S_OLE1_FAR_RM
#myristic_14_0_percent        0.08697556    0.6586564      0.3207597               0.40597013
#palmitoleic_16_1_percent     0.60211756    0.2662922      0.2617028               0.04404311
#palmitic_16_0_percent        0.05719410    0.6354966      0.1801385               0.59184509
#oleic_18_1_percent           0.01126871    0.7430847      0.2584599               0.39601910
#stearic_18_0_percent         0.35806366    0.2266299      0.4239089               0.21219673
#totalSaturated_percent       0.11832352    0.4191694      0.2542052               0.44281148
#totalUnsaturated_percent     0.12177791    0.4158889      0.2510617               0.37797694
#total16_percent              0.003526482   0.004443969    0.2822321               0.10539990
#total18_percent              0.004808213   0.005670757    0.1968508               0.10926480

# estimates
#                                    RM           BY_OAF1_L63S   BY_OLE1_FAR_RM          BY_OAF1_L63S_OLE1_FAR_RM
#myristic_14_0_percent.mean of y     0.004267141  0.001759857    0.004255007             -0.002543336
#palmitoleic_16_1_percent.mean of y  0.020886383 -0.076406825   -0.051083133             -0.077415402
#palmitic_16_0_percent.mean of y    -0.082463979  0.029403944    0.065894629              0.018328608
#oleic_18_1_percent.mean of y        0.075438034 -0.014848193   -0.047522073              0.033412771
#stearic_18_0_percent.mean of y     -0.018870741  0.059327079    0.028752835              0.023066865
#totalSaturated_percent.mean of y   -0.097067579  0.090490880    0.098902471              0.038852137
#totalUnsaturated_percent.mean of y  0.096324416 -0.091255018   -0.098605206             -0.044002631
#total16_percent.mean of y          -0.061577596 -0.047002881    0.014811496             -0.059086794
#total18_percent.mean of y           0.056567293  0.044478886   -0.018769238              0.056479636



pdf("allCombined_composition.pdf", width=11, height=6)
print(ggplot(datGather[datGather$genotype !="RM" & !str_detect(datGather$compound, "total"),], aes(x=compound, y=fraction, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

print(ggplot(datGather[!str_detect(datGather$compound, "total"),], aes(x=compound, y=fraction, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75)) +
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF", "purple")) + theme_light())

print(ggplot(datGather[datGather$genotype !="RM" & str_detect(datGather$compound, "total"),], aes(x=compound, y=fraction, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

print(ggplot(datGather[str_detect(datGather$compound, "total"),], aes(x=compound, y=fraction, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF", "purple")) + theme_light())

print(ggplot(datGather[datGather$genotype !="RM" & str_detect(datGather$compound, "totalUnsaturated"),], aes(x=compound, y=fraction, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

print(ggplot(datGather[datGather$genotype !="RM" & str_detect(datGather$compound, "total18"),], aes(x=compound, y=fraction, color=genotype)) +
geom_jitter(position=position_jitterdodge()) +
stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1), geom="pointrange", position=position_dodge(0.75))+
scale_color_manual(values=c("grey", "salmon", "#66B2FF", "#CC99FF")) + theme_light())

dev.off()
