
CUT SCRIPTS


```{r, eval = FALSE}
# - source AA glycine, serine, and threonine removed, showing 13C.CSIA.no.glyserthr


##################
# source AA glycine and serine removed also threonine
dfC.trim<-d13C.dat.long[!(d13C.dat.long$AA.short=="Ser" | 
                            d13C.dat.long$AA.short=="Gly" | 
                            d13C.dat.long$AA.short=="Thr"),]

AA.means<-aggregate(d13C.value~AA.cat+Treat.Int+Fraction, data=dfC.trim, mean, na.rm=TRUE)
AA.sd<-aggregate(d13C.value~AA.cat+Treat.Int+Fraction, data=dfC.trim, na.rm=TRUE, sd)
colnames(AA.sd)[4]="SD"
AA.means<-cbind(AA.means, AA.sd[4])
AA.means$Fraction<-factor(AA.means$Fraction, levels=c("host", "symb", "plank"))

pd <- position_dodge(0.5) #offset for error bars
ggplot(AA.means, aes(x=Treat.Int, y=d13C.value, color=Fraction, shape=AA.cat)) +
  geom_errorbar(aes(ymin=d13C.value-SD, ymax=d13C.value+SD, color=Fraction, shape=AA.cat), size=.5, width=0, position=pd) +
  geom_point(aes(color=Fraction, shape=AA.cat), size=3, position=pd) +
  coord_cartesian(ylim=c(-30, -5)) + 
  scale_x_discrete(name ="Treatments", 
                   labels=c("Dark\nFed", "Light\nFed", "Light\nNot Fed", "Plankton")) +
  ggtitle(expression(paste(delta^{13}, C[AA], " by fraction, Tr/So-AA, - Ser/Gly/Thr"))) +
  scale_color_manual(values=c("coral", "springgreen3", "skyblue3")) +
  ylab(expression(paste(delta^{13}, C[AA], " (\u2030, V-PDB)"))) +
  Fig.formatting
ggsave("figures/carbon/d13C.CSIA.no.glyserthr.pdf", height=5, width=8, encod="MacRoman")
```


```{r, eval=FALSE}
#threonine in host and symbiont (microbial source), showing d13C.CSIA_TrSo.thr
#

#################
# threonine in host and symbiont (microbial source)
dfC.thr<-d13C.dat.long[(d13C.dat.long$AA.short=="Thr"),]

AA.means<-aggregate(d13C.value~AA.short+Treat.Int+Fraction, data=dfC.thr, mean, na.rm=TRUE)
AA.sd<-aggregate(d13C.value~AA.short+Treat.Int+Fraction, data=dfC.trim, na.rm=TRUE, sd)
colnames(AA.sd)[4]="SD"
AA.means<-cbind(AA.means, AA.sd[4])
AA.means$Fraction<-factor(AA.means$Fraction, levels=c("host", "symb", "plank"))

ggplot(AA.means, aes(x=Treat.Int, y=d13C.value, color=Fraction)) +
  geom_errorbar(aes(ymin=d13C.value-SD, ymax=d13C.value+SD), size=.5, width=0, position=pd) +
  geom_point(aes(color=Fraction), size=3, position=pd) +
  coord_cartesian(ylim=c(0, -30)) + 
  scale_x_discrete(name ="Treatments", 
                   labels=c("Dark\nFed", "Light\nFed", "Light\nNot Fed", "Plankton")) +
  ggtitle(expression(paste("Threonine ", delta^{13}, C[AA]))) +
  scale_color_manual(values=c("coral", "springgreen3", "skyblue3")) +
  ylab(expression(paste(delta^{13}, C[AA], " (\u2030, v-PDB)"))) +
  Fig.formatting
ggsave("figures/carbon/d13C.CSIA_TrSo.thr.pdf", height=5, width=8, encod="MacRoman")
```


```{r, eval =FALSE}
# all AA pooled by fraction
pd <- position_dodge(0.5) #offset for error bars
ggplot(df.mean, aes(x=AA.short, y=d13C.value)) +
  
  
  geom_errorbar(aes(ymin=d13C.value-SD, ymax=d13C.value+SD, color=Fraction), 
                size=.5, width=0, position=pd) +
  geom_point(size=3, pch=19,  position=pd, aes(color=Fraction)) +
  geom_vline(xintercept=7.5, linetype="solid", color = "gray") +
  annotate(geom="text", label="Trophic-AA", x=4.0, y=0, color="gray40") +
  annotate(geom="text", label="Source-AA", x=10.0, y=0, color="gray40") +
  ggtitle(expression(paste(delta^{13}, C[AA], " by biological fraction"))) +
  coord_cartesian(ylim=c(-35, 0)) + 
  xlab("Amino Acids") + 
  ylab(expression(paste(delta^{13}, C[AA], " (\u2030, V-PDB)"))) +
  scale_color_manual(values=c("coral", "springgreen3", "skyblue3")) +
  Fig.formatting
ggsave("figures/carbon/d13C.CSIA_Fraction.pdf", height=5, width=8, encod="MacRoman")
```


df.mean<-aggregate(d13C.value~AA.short+Fraction, data=dfC, mean, na.rm=TRUE)
df.n<-aggregate(d13C.value~AA.short+Fraction, data=dfC, length)
df.SD<-aggregate(d13C.value~AA.short+Fraction, data=dfC, sd, na.rm=TRUE)
colnames(df.SD)[3]="SD"
df.mean<-cbind(df.mean, df.SD[3])
df.mean$Fraction<-factor(df.mean$Fraction, levels=c("host", "symb", "plank"))


- threonine in host and symbiont (microbial source), showing d15N.CSIA_TrSo.thr
```{r, results='hide', warning=FALSE}
##################
# threonine in host and symbiont (microbial source)
dfN.thr<-d15N.dat.long[(d15N.dat.long$AA.short=="Thr"),]
AA.means<-aggregate(d15N.value~AA.cat+AA.short+Treat.Int+Fraction, data=dfN.thr, mean, na.rm=TRUE); AA.means
AA.sd<-aggregate(d15N.value~AA.cat+Treat.Int+Fraction, data=dfN.trim, na.rm=TRUE, sd)
colnames(AA.sd)[4]="SD"
AA.means<-cbind(AA.means, AA.sd[4])
AA.means$Fraction<-factor(AA.means$Fraction, levels=c("host", "symb", "plank"))

ggplot(AA.means, aes(x=Treat.Int, y=d15N.value, color=Fraction)) +
  geom_errorbar(aes(ymin=d15N.value-SD, ymax=d15N.value+SD), size=.5, width=0, position=pd) +
  geom_point(aes(color=Fraction), size=3, position=pd) +
  coord_cartesian(ylim=c(6, -6)) + 
  scale_x_discrete(name ="Treatments", 
                   labels=c("Dark\nFed", "Light\nFed", "Light\nNot Fed", "Plankton")) +
  ggtitle(expression(paste("Threonine ", delta^{15}, N[AA]))) +
  scale_color_manual(values=c("coral", "springgreen3", "skyblue3")) +
  ylab(expression(paste(delta^{15}, N[AA], " (\u2030, Air)"))) +
  Fig.formatting
ggsave("figures/nitrogen/d15N.CSIA_TrSo.thr.pdf", height=5, width=8, encod="MacRoman")
```



- source and trophic AA w/o glycine, serine and threonine, showing d15N.CSIA_TrSo.no.glyserthr
```{r, results='hide', warning=FALSE}
##################
# source and trophic AA w/o glycine  serine and threonine
dfN.trim<-d15N.dat.long[!(d15N.dat.long$AA.short=="Ser" | 
                            d15N.dat.long$AA.short=="Gly" | d15N.dat.long$AA.short=="Thr"), ]

AA.means<-aggregate(d15N.value~AA.cat+Treat.Int+Fraction, data=dfN.trim, mean, na.rm=TRUE)
AA.sd<-aggregate(d15N.value~AA.cat+Treat.Int+Fraction, data=dfN.trim, na.rm=TRUE, sd)
colnames(AA.sd)[4]="SD"
AA.means<-cbind(AA.means, AA.sd[4])
AA.means$Fraction<-factor(AA.means$Fraction, levels=c("host", "symb", "plank"))


ggplot(AA.means, aes(x=Treat.Int, y=d15N.value, shape=AA.cat, color=Fraction)) +
  geom_errorbar(aes(ymin=d15N.value-SD, ymax=d15N.value+SD, color=Fraction, shape=AA.cat), size=.5, width=0, position=pd) +
  geom_point(aes(color=Fraction, shape=AA.cat), size=3, position=pd) +
  coord_cartesian(ylim=c(14, 0)) + 
  scale_x_discrete(name ="Treatments", 
                   labels=c("Dark\nFed", "Light\nFed", "Light\nNot Fed")) +
  ggtitle(expression(paste(delta^{15}, N[AA], " by biological fraction, Tr/So-AA, - Ser/Gly/Thr"))) +
  scale_color_manual(values=c("coral", "springgreen3",  "skyblue3")) +
  ylab(expression(paste(delta^{15}, N[AA], " (\u2030, Air)"))) +
  Fig.formatting
ggsave("figures/nitrogen/d15N.CSIA_TrSo.no.glyserthr.pdf", height=5, width=8, encod="MacRoman")
```

