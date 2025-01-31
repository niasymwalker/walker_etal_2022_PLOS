#R Code used in the analysis of Walker et al. 2022 PLOS
#Titled: "Persistence of phenotypic responses to short-term heat stress in the tabletop coral Acropora hyacinthus"
#Contact Nia Walker at niawalker13@gmail.com for any further inquiries.

#Packages used in the analysis:
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(nlme)
library(lme4)
library(sjstats)
library(png)
library(cowplot)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(vcrpart)
library(grid)

#Datasets
Maindata <- read_excel("./data/maindata.xlsx", na = "NA")
Symprops <- read_excel("./data/symprops.xlsx", na = "NA")
CoralVBSphoto <- readPNG("./data/CoralVBSphoto.png")
ReefRegions <- read_excel("./data/reefregions.xlsx", na = "NA")

##################Figure 1 Analysis: Symbiont Density vs. VBS Scores#####
Symprops$VBSValues <- as.factor(Symprops$VBSValues)

Fig1 <- ggplot(Symprops, aes(x = VBSValues, y = SymProp)) +
  geom_point(aes(x=VBSValues,y=SymProp,shape=ExpGroups,color=ExpGroups),data=Symprops, position=position_jitter(w=0.25,h=0), size=3) +
  geom_boxplot(aes(x=VBSValues,y=SymProp), data=Symprops, alpha=0) +
  scale_color_manual(values=c("dodgerblue4", "green4")) + 
  annotate("text", x=0.7, y=0.095, label="a", size=4) +
  annotate("text", x=1.7, y=0.126, label="a", size=4) +
  annotate("text", x=2.7, y=0.078, label="b", size=4) +
  annotate("text", x=3.65, y=0.039, label="c", size=4) +
  annotate("text", x=4.68, y=0.046, label="c", size=4) +
  labs(color="Heat Stress\nExperiment\nCategories",shape="Heat Stress\nExperiment\nCategories", x = "Visual Bleaching Score (VBS)",
       y = "Symbiont Concentration\n(Symbiont Cell Counts/Total Cell Counts)") + 
  theme_bw() + theme(axis.title.y = element_text(margin = margin(l = 10, r=5),hjust=0.5,vjust=0.5),
                         axis.text.y = element_text(hjust=0.5,vjust=0.5),
                         legend.text = element_text(
                           hjust = 0.5,         # Center the text horizontally
                           vjust = 0.5          # Center the text vertically
                         ))

image_plot <- ggdraw() +
  draw_image(CoralVBSphoto, scale=1)

quartz(w=7.2,h=6)
plot_grid(Fig1,image_plot,ncol=1,align="v",rel_heights = c(2.5,1))
quartz.save("./output/Fig1.pdf",type="pdf")

#tukey test for all VBS Scores 1-5, plus random region variable

guavagroups_nlme <- lme(SymProp ~ VBSValues, random=~1|Region,method="ML",data=Symprops)
guavagroups_anova2 <- anova(lme(SymProp ~ VBSValues, random=~1|Region, method = "ML", data=Symprops))
summary(glht(guavagroups_nlme, linfct = mcp(VBSValues = "Tukey")))
r2(guavagroups_nlme)

#tukey test for VBS groups: VBS < 3, VBS = 3, VBS > 3
Symprops$ResistanceGroups <- as.factor(Symprops$ResistanceGroups)
guavaresist_nlme <- lme(SymProp ~ ResistanceGroups, random=~1|Region,method="ML",data=Symprops)
guavaresist_anova2 <- anova(lme(SymProp ~ ResistanceGroups, random=~1|Region, method = "ML", data = Symprops))
summary(glht(guavaresist_nlme, linfct = mcp(ResistanceGroups = "Tukey")))
r2(guavaresist_nlme)

##################Figure 2 Analysis#####
#Convert data from wide to long, ALL WIDE COLUMNS (VBS, Mort, Algae)
Maindata.long <- Maindata %>%
  gather(VBSDay, VBS, VBSDay_0:VBSDay_7) %>%
  gather(CtrlVBSDay, CtrlVBS, CtrlVBSDay_0:CtrlVBSDay_7) %>%
  gather(MortDay, Mort, MortDay_0:MortDay_7) %>%
  gather(VBSDiffDay, VBSDiff, VBSDiffDay_0:VBSDiffDay_7) %>%
  gather(CtrlMortDay, CtrlMort, CtrlMortDay_0:CtrlMortDay_7) %>%
  distinct()
#convert to factor labels (VBS, Mort, Algae)
#Daily.FINAL <- mutate_at(Daily.FINAL, vars('VBS', "Mort", "CtrlVBS", "CtrlMort"), as.factor())
Factor_VBS <- factor(Maindata.long$VBS, labels = c("none", "none-visible", "visible", "visible-moderate", "moderate", "moderate-severe", "severe", "total"))
Factor_CtrlVBS <- factor(Maindata.long$CtrlVBS, labels = c("none", "none-visible", "visible", "visible-moderate", "moderate"))
Factor_Mort <- factor(Maindata.long$Mort, labels = c("D", "PD", "ND"))
Factor_CtrlMort <- factor(Maindata.long$CtrlMort, labels = c("D", "PD", "ND"))
Factor_Coral <- factor(Maindata.long$Coral)
Factor_Resistance <- factor(Maindata.long$Resistance, levels = c("Low", "Moderate", "High"), labels = c("Low Bleaching", "Moderate Bleaching", "High Bleaching"))
#Factor_Day <- factor(Daily.FINAL$Day, labels = c("0", "1", "2", "3", "4", "5", "6", "7"))
Factor_Diff <- factor(Maindata.long$VBSDiff, labels = c("-1.5", "1", "-0.5", "0", "0.5", "1", "1.5", "2", "2.5", "3", "4"))
Factor_VBSDay <- factor(Maindata.long$VBSDay, labels = c("0", "1", "2", "3", "4", "5", "6", "7"))
Factor_CtrlVBSDay <- factor(Maindata.long$CtrlVBSDay, labels = c("0", "1", "2", "3", "4", "5", "6", "7"))
Factor_MortDay <- factor(Maindata.long$MortDay, labels = c("0", "1", "2", "3", "4", "5", "6", "7"))
Factor_VBSDiffDay <- factor(Maindata.long$VBSDiffDay, labels = c("0", "1", "2", "3", "4", "5", "6", "7"))
Factor_CtrlMortDay <- factor(Maindata.long$CtrlMortDay, labels = c("0", "1", "2", "3", "4", "5", "6", "7"))
Factor_VBSDay0 <- factor(Maindata$VBSDay_0, labels = c("none", "none-visible", "visible", "visible-moderate", "moderate", "moderate-severe", "severe", "total"))

FactOrd_VBS <- as.ordered(Factor_VBS)
FactOrd_CtrlVBS <- as.ordered(Factor_CtrlVBS)
FactOrd_Mort <- as.ordered(Factor_Mort)
FactOrd_CtrlMort <- as.ordered(Factor_CtrlMort)
FactOrd_Diff <- as.ordered(Factor_Diff)
#FactOrd_Day <- as.ordered(Factor_Day)
FactOrd_VBSDay <- as.ordered(Factor_VBSDay)
FactOrd_CtrlVBSDay <- as.ordered(Factor_CtrlVBSDay)
FactOrd_MortDay <- as.ordered(Factor_MortDay)
FactOrd_CtrlMortDay <- as.ordered(Factor_CtrlMortDay)
FactorOrd_VBSDiffDay <- as.ordered(Factor_VBSDiffDay)
FactOrd_VBSDay0 <- as.ordered(Factor_VBSDay0)
#dataframe for factor variables

FACTORED_Maindata.long <- data.frame(Maindata.long$Reef, Maindata.long$VBSDay, Maindata.long$VBSDiffDay, Maindata.long$CtrlVBSDay, Maindata.long$MortDay, Maindata.long$CtrlMortDay, Factor_Coral, 
                                     Factor_Resistance, FactOrd_VBS, FactOrd_CtrlVBS, FactOrd_Mort, FactOrd_CtrlMort, FactOrd_VBSDay, FactOrd_CtrlVBSDay, FactOrd_MortDay, FactOrd_CtrlMortDay, 
                                     FactorOrd_VBSDiffDay, FactOrd_VBSDay0)

#bleaching score ordinal timeseries trajectory
par(mfrow=c(3,1))

quartz(w=7,h=6)
otsplot(x = FACTORED_Maindata.long$FactOrd_VBSDay, y = FACTORED_Maindata.long$FactOrd_VBS, 
        subject = FACTORED_Maindata.long$Factor_Coral, groups = FACTORED_Maindata.long$Factor_Resistance, 
        filter = otsplot_filter("cumfreq", level = 0.5), cex = 1, lwd=10, xlab="Days Post-Stress")

quartz.save("./output/Fig2.pdf",type="pdf")
###NOTE: edited cutoff y-axis, enlarged the lines, and added rectangles manually; I was not able to do this in R.

##################Figure 3 Analysis#####
#Figure 3a: VBS on Day 0 vs. Day 7
Maindata$Resistance <- factor(Maindata$Resistance, levels = c("High", "Moderate", "Low"))
M3A <- ggplot(Maindata, aes(x=VBSDay_0, y=VBSDay_7)) +
  geom_point(data=Maindata, aes(group=Resistance,color=Resistance,shape=Resistance, fill=Resistance), size=4, position=position_jitter(w=.05,h=.05), color = "black", alpha=.8) +
  scale_shape_manual(values=c(21,24,22)) + 
  scale_fill_manual(values=c("coral3", "goldenrod1", "skyblue3")) +
  scale_color_manual(values=c("coral3", "goldenrod1", "skyblue3")) +
  geom_smooth(method = "lm",
              se = TRUE, level = 0.95) +
  labs(x = "VBS on Day 0 of Recovery",
       y = expression("VBS on Day 7 of Recovery")) + theme_bw() + theme(axis.text = element_text(size=10), axis.title = element_text(size=10), legend.text = element_text(size=10), legend.title = element_text(size=10)) 

VBSDay0 <- c(1,2,1.5,2,2,2,2,2,2,2,2,2.5,2,2,2.5,2.5,2.5,2.5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3.5,4,3,4,4,4,4,4,4,4,4,4)
VBSDay7 <- c(1,1,1.5,2,2,2,2,2,2,2,2,2,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3.5,3.5,3.5,4,4,4,4,4,4,4)
ReefRegion <- c("SL", "NL", "FR", "SL","SL","WSL","WSL","WSL","WSL","NL","NL","WSL","EL","EL","FR","NL","FR","FR","FR","FR","EL","EL","EL","FR","FR","FR","SL","SL","SL","SL","WSL","WSL","NL","NL","NL","NL","NL","NL","SL","NL","EL","NL","NL","EL","FR","FR","NL","NL","NL","NL")
VBSdataframe <- data.frame(VBS0,VBS7,ReefRegion)

VBS0v7_nlme <- lme(VBSDay7 ~ VBSDay0, random=~1|ReefRegion,method="ML", data=Dataset)
r2(VBS0v7_nlme)

#Figure 3b: Mortality from Recovery Day 0-7, separated by bleaching severity group
#mortality bar plot done as percentages per group
Severity <- c("High", "Moderate", "Low", "High", "Moderate", "Low", "High", "Moderate", "Low", "High", "Moderate", "Low", "High", "Moderate", "Low", "High", "Moderate", "Low", "High", "Moderate", "Low", "High", "Moderate", "Low","High", "Moderate", "Low","High", "Moderate", "Low")
Severity <- factor(Severity, levels = c("High", "Moderate", "Low"))

MDay <- c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,10,10,10,11,11,11)
Mortality <- c(31.71,0,0,46.34,7.89,0,60.98,31.58,13.64,68.29,39.47,13.64,73.17,42.11,13.64,73.17,42.11,13.64,73.17,42.11,18.18,73.17,44.74,18.18,73.17,44.74,18.18,73.17,44.74,22.72)

barmortality_corrected <- data.frame(Severity, MDay, Mortality)

M3B <- ggplot(barmortality_corrected, aes(x=MDay, fill=Severity, y=Mortality)) +
  geom_bar(position="dodge",stat="identity",color="black") +
  scale_color_manual(values=c("red4", "goldenrod4", "blue4")) +
  scale_fill_manual(values=c("coral3", "goldenrod1", "skyblue3")) +
  geom_vline(xintercept = 2.5, linetype = "solid", color="black") +  
  geom_vline(xintercept = 8.5, linetype = "dotted", color="black") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,10,11), labels=c("0","1","2","3","4","5","6","7","~30","~60")) +
  labs(fill="Bleaching Severity\nCategories", x="Days Post-Heat Stress (Recovery Period)", y=expression("Bleached Coral Sample Mortality (%)")) +
  ylim(0,100) + theme_bw() + theme(axis.text = element_text(size=10), axis.title = element_text(size=10), legend.text = element_text(size=10), legend.title = element_text(size=10))


#Figure 3c Analysis: VBS day 0 vs. Mort Day 60
#VBS Day 0 vs. Mortality Day 60
VBSDay_0_Logit <- c(3,4,3,3,3,2,2,3,4,3,2,3,3.5,3,4,3,3,4,3,3,3,4,3,3,2,1,3,2,3,3,3,3,2,4,2,2,3,4,2.5,2,2,3,3,3,4,3,3,2.5,1.5,4,4,2,3,4,4,3,3,4,3,2,4,5,4,2,4,4,4,4,2,3,2.5,3,3,3,5,4,4,3,5,4,5,5,3,3,2.5,2,3,2.5)
VBSDay_0_Resist <- c("Moderate","High","Moderate","Moderate","Moderate","Low","Low","Moderate","High","Moderate","Low","Moderate","High","Moderate","High","Moderate","Moderate","High","Moderate","Moderate","Moderate","High","Moderate","Moderate","Low","Low","Moderate","Low","Moderate","Moderate","Moderate","Moderate","Low","High","Low","Low","Moderate","High","Low","Low","Low","Moderate","Moderate","Moderate","High","Moderate","Moderate","Low","Low","High","High","Low","Moderate","High",
                     "High","Moderate","Moderate","High","Moderate","Low","High","High","High","Low","High","High","High","High","Low","Moderate","Low","Moderate","Moderate","Moderate","High","High","High","Moderate","High","High","High","High","Moderate","Moderate","Low","Low","Moderate","Low")
AliveDay_60_Logit <- c(0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1,0,1,1,0,0,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,0,1,1,0,0,0,1,1,0)
ReefRegion_Logit <- c("EL","EL","EL","EL","EL","EL","EL","EL","EL","EL","EL","SL","SL","FR","FR","FR","FR","FR","FR","FR","FR","SL","SL","SL","SL","SL","SL","SL","SL","SL","SL","SL","SL","SL","WSL","WSL","WSL","WSL","WSL","WSL","WSL","WSL","WSL","FR","FR","FR","FR","FR","FR","FR","FR","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","NL","FR","FR","FR","FR","FR")

df_logitVBSMort <- data.frame(VBSDay_0_Logit, AliveDay_60_Logit,ReefRegion_Logit, VBSDay_0_Resist)
logit_VBSMort60 <- glmer(AliveDay_60_Logit ~ VBSDay_0_Logit + (1|ReefRegion_Logit), data=df_logitVBSMort, family = binomial)
summary(logit_VBSMort60)
r2(logit_VBSMort60)

M3C <- ggplot(df_logitVBSMort, aes(x = VBSDay_0_Logit, y = AliveDay_60_Logit,color=VBSDay_0_Resist)) +
  geom_count(aes(size = after_stat(prop), group = VBSDay_0_Logit), shape = 21, fill = c("white"), stroke = 2) +
  #scale_size_area(max_size=20) +
  scale_color_manual(values=c("coral3", "goldenrod1", "skyblue3")) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE, color="black") +
  labs(x = "Visual Bleaching Score (VBS)\nDay 0 of Recovery",
       y = "Probability of Total Mortality\n2 Months of Recovery", size ="Proportion of Living\nvs. Dead Samples by\nBleaching Severity") + 
  guides(color="none") +
  xlim(1,5) + ylim(0,1) + theme_bw() + theme(axis.text = element_text(size=10), axis.title = element_text(size=10), legend.text = element_text(size=10), legend.title = element_text(size=10))

quartz(w=10,h=6.5)

plot_grid(M3A,M3B,M3C,NULL, labels=c("A","B","C",""), align="hv", axis="tb", ncol=2, rel_widths = c(1,1,1,1))
quartz.save("./output/Fig3.pdf",type="pdf")

##################Figure 4 Analysis#####
#Figure 4a: Bleaching severity across variable reef environments, stacked barplot 
ReefRegions$Reef_Prop <- as.factor(ReefRegions$Reef_Prop)
ReefRegions$Prop_Group <- factor(ReefRegions$Prop_Group, levels = c("Low", "Moderate","High"))

M4A <- ggplot(ReefRegions, aes(x=Reef_Prop, y= Props, fill = Prop_Group, color = "black")) +
  geom_bar(position="stack", stat="identity", color = "black", width=0.7) +
  scale_fill_manual(values=c("skyblue3","goldenrod1", "coral3")) + 
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  annotate("text", x=c(5.75),y=c(1.045), label=c("Southern\nPatch Reefs"),size=3) +
  geom_vline(xintercept = 11.5, linetype="solid",color="black")+
  annotate("text", x=c(16.5),y=c(1.045), label=c("Northern\nPatch Reefs"),size=3) +
  geom_vline(xintercept = 21.5, linetype="solid",color="black")+
  annotate("text", x=c(23.5),y=c(1.045), label=c("South\nFRs"),size=3) + 
  geom_vline(xintercept = 25.5, linetype="solid",color="black")+
  annotate("text", x=c(27),y=c(1.045), label=c("North\nFRs"),size=3) +
  theme(aspect.ratio = 7/3) +
  labs(x="Individual Reef Sites", y="Proportion of Colonies in Severity Categories") +
  theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

#Figure 4b: Bleaching severity across variable reef environments, boxplot
BleachingSeverity_Cats <- c("High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Dead","Dead", "Dead","Dead", "Dead","Dead", "Dead","Dead")
BleachingSeverity_Cats <- factor(BleachingSeverity_Cats, levels = c("Low", "Moderate", "High", "Dead"))
Reef_Temps <- c(1845,1547,2451.666667,907.5,907.5,1425,1412.666667,3101.75,1221,1092.5,1918,495.5,1127.5,1127.5,1933.5,871,871,871,1845,1845,1845,1547,1547,2451.666667,907.5,907.5,907.5,907.5,907.5,907.5,2652.333333,2652.333333,2652.333333,2652.333333,1204,3101.75,2395,2395,1221,1221,1221,495.5,1933.5,871,151.333333,1515,1547,1547,2652.333333,1412.666667,3101.75,3101.75,3101.75,2395,2395,495.5,1526,1526,1515,1425,1204,1221,1221,1221,1092.5,1092.5,1918)
tempvsbleaching <- data.frame(BleachingSeverity_Cats,Reef_Temps)  
#tempvsbleaching$BleachingSeverity_Cats <- as.factor(tempvsbleaching$BleachingSeverity_Cats)

M4B <- ggplot(tempvsbleaching, aes(x=BleachingSeverity_Cats, y=Reef_Temps, fill=BleachingSeverity_Cats)) +
  geom_boxplot() +
  geom_jitter(width=0.1, alpha=0.8, size=1, color="black") +
  scale_fill_manual(values=c("skyblue3", "goldenrod1", "coral3", "gray40")) +
  labs(fill="Bleaching Severity\nCategories",x="Bleaching Severity Categories", y=expression("Temperature (10 Min Interval Counts Above 31"*degree*C*")")) +
  theme_bw() + theme(legend.position = "left")

quartz(w=10,h=5)
plot_grid(M4A,M4B, labels=c("A","B"), ncol=2, rel_widths = c(1,1.5), label_x = c(0, 0.25))
quartz.save("./output/Fig4.pdf", type="pdf")

tempbleaching <- polr(BleachingSeverity_Cats ~ Reef_Temps, data=tempvsbleaching, Hess=TRUE)
ctable <- coef(summary(tempbleaching))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
r2(tempbleaching)

##################Figure 5 Analysis#####
#keep figure as is for Fig 5 day0,7,30,60
vbs_projections <- c(3.5,3,3,2.5,3,4,2.5,3,3,4,3,2.5,3,3,2,1,3,2,3,3,2,2,2,2,2,3,3,2.5,1.5,4,3,3,3,4,2,3.5,1,4,3,4,4,2,2.5,3,3,3.5,3,2.5,2.5,2.5,2,1,1,1,2.5,2,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,2,2,1,1,1,1,1,2.5,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,NA,1,1,1,1,1,1,1,1)
resist_projections <- c("Moderate","Moderate","Moderate","Low","Moderate","High","Low","Moderate","Moderate","High","Moderate","Moderate","Moderate","High","Low","Low","Moderate","Low","Moderate","Moderate","Low","Low","Low","Low","Low","Moderate","Moderate","Low","Low","High","Moderate","Moderate","Moderate","High","Low","High","Low","High","High","High","High","Low","Low","Moderate","Moderate","High","Moderate","Moderate","Low","Low","Moderate","Moderate","Moderate","Low","Moderate","High","Low","Moderate","Moderate","High","Moderate","Moderate","Moderate","High","Low","Low","Moderate","Low","Moderate","Moderate","Low","Low","Low","Low","Low","Moderate","Moderate","Low","Low","High","Moderate","Moderate","Moderate","High","Low","High","Low","High","High","High","High","Low","Low","Moderate","Moderate","High","Moderate","Moderate","Low","Low","Moderate","Moderate","Moderate","Low","Moderate","High","Low","Moderate","Moderate","High","Moderate","Moderate","Moderate","High","Low","Low","Moderate","Low","Moderate","Moderate","Low","Low","Low","Low","Low","Moderate","Moderate","Low","Low","High","Moderate","Moderate","Moderate","High","Low","High","Low","High","High","High","High","Low","Low","Moderate","Moderate","High","Moderate","Moderate","Low","Low")
resist_projections <- factor(resist_projections, levels = c("High","Moderate","Low"))
day_projections <- c("Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 7","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 30","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60","Day 60")
day_projections <- factor(day_projections, levels = c("Day 7","Day 30","Day 60"))
mort_projections <- c("ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","PD","ND","PD","PD","ND","PD","PD","PD","ND","ND","ND","ND","ND","PD","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","PD","ND","ND","ND","ND","ND","PD","ND","ND","PD","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","PD","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","ND","PD","PD","ND","ND","ND","ND","ND","ND","PD","D","ND","PD","ND","ND","ND","ND","ND","ND")
coral_projections <- c(3,6,7,13,18,20,25,35,41,46,47,48,50,66,69,70,71,72,74,75,101,104,110,135,137,138,139,171,179,180,207,224,227,232,236,251,252,255,257,258,263,266,268,270,273,332,340,373,383,394,3,6,7,13,18,20,25,35,41,46,47,48,50,66,69,70,71,72,74,75,101,104,110,135,137,138,139,171,179,180,207,224,227,232,236,251,252,255,257,258,263,266,268,270,273,332,340,373,383,394,3,6,7,13,18,20,25,35,41,46,47,48,50,66,69,70,71,72,74,75,101,104,110,135,137,138,139,171,179,180,207,224,227,232,236,251,252,255,257,258,263,266,268,270,273,332,340,373,383,394)
coral_projections <- as.factor(coral_projections)
projections_df <- data.frame(vbs_projections,resist_projections,day_projections,mort_projections,coral_projections, na.rm = TRUE)

M5 <- ggplot(projections_df, aes(x = day_projections, y = vbs_projections)) +
  geom_point(aes(x=day_projections,y=vbs_projections,shape=resist_projections,color=resist_projections, group=resist_projections),data=projections_df, position=position_jitterdodge(jitter.width=0.5,jitter.height=0.1,dodge.width=0.75), size=1.5) +
  geom_violin(aes(color=resist_projections, shape=resist_projections), scale = "count", alpha=0) +
  geom_vline(xintercept = c(1.5), linetype="solid", color="black") +
  geom_vline(xintercept = c(2.5), linetype="solid", color="black") +
  scale_color_manual(values=c("coral3", "goldenrod1", "skyblue3")) +
  annotate("text", x=0.7, y=4.2, label="a", size=4) +
  annotate("text", x=1, y=3.7, label="b", size=4) +
  annotate("text", x=1.3, y=2.7, label="c", size=4) +
  labs(shape = "Bleaching Severity\nCategories", color = "Bleaching Severity\nCategories", x = "Days Post-Heat Stress (Recovery Period)",
       y = expression("Visual Bleaching Score (VBS)")
  ) + theme_bw()

quartz(w=7,h=5)
plot_grid(M5)
quartz.save("./output/Fig5.pdf", type="pdf")
