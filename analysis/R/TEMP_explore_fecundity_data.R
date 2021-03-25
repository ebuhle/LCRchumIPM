#Temp explore fecundity data
library(tidyverse)
library(ggplot2)
dat<-read_csv("data/Data_ChumFecundity_fromHatcheryPrograms_2017-01-25.csv")%>%
  mutate(EggMass =`Green egg avg weight` * `Estimated Fecundity`,
         BodyWeight = EggMass/(`Reproductive Effort` / 100),
         PartialSpawn = factor(ifelse(`Reproductive Effort`<16,1,0))
         )%>%
  filter(Age%in%c(3:5),
         BodyWeight <10000
         )

ggplot(dat,aes(y=`Estimated Fecundity`,x=BodyWeight,color=PartialSpawn))+
  geom_point(alpha=0.3)+
  facet_wrap(~Age)

quants<-c(0.025, 0.25, 0.5, 0.75, 0.975)
dat%>%
  group_by(Age,PartialSpawn)%>%
  summarise(Fecundity = quantile(`Estimated Fecundity`, quants), q = quants)%>%
  as.data.frame()%>%
  pivot_wider(names_from=q,values_from=Fecundity)
  print()
