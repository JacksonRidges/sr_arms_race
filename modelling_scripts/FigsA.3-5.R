#Load libraries
library(ggplot2)
library(cowplot)

#Define parameter values (see README for definitions)
d=1/2; dsup=0; ssr=0.59; h=0.135; ssup=0.2; ssry=0.2; ssrsup=ssry;

#Create a blank data frame for results
df=data.frame(dsup=numeric(),ssup=numeric(),ssry=numeric(),xsreqm=numeric(),ysup=numeric())

#Calculate equilibrium frequencies (see Mathematica file)
for(dsup in c(0)) { #assume no drive if suppressed
  for(ssup in seq(0.01,0.6,0.005)) { #range of costs of suppression
    for(ssry in seq(0.01,0.5,0.005)) { #range of costs of drive in males
        ssrsup=1-(1-ssup)*(1-ssry)
        df[length(df$dsup)+1,]=c(dsup,ssup,ssry,
          1-((2*dsup*(-1+ssrsup)-ssrsup-2*d*(-1+ssry)+
          ssry)*(-ssup-(-2+h)*ssup*ssr+2*dsup*(-1+h*ssr)*(-1+ssrsup)+
          ssrsup-2*d*(-1+h*ssr)*(-1+ssry)-ssry+
          h*ssr*(-ssrsup+ssry)))/(-ssup^2+h*ssup^2*ssr+
          4*dsup^2*(-1+h*ssr)*(-1+ssrsup)^2+2*ssup*ssrsup-2*ssup*ssr*ssrsup+
          2*h*ssup*ssr*ssrsup-ssrsup^2+h*ssr*ssrsup^2-4*dsup*(-1+
          ssrsup)*(ssup+(-1+h)*ssup*ssr+(-1+h*ssr)*(ssrsup-ssry)+
          2*d*(-1+h*ssr)*(-1+ssry))+
          4*d*(ssup+(-1+h)*ssup*ssr+(-1+h*ssr)*(ssrsup-ssry))*(-1+
          ssry)+4*d^2*(-1+h*ssr)*(-1+ssry)^2-
          2*(ssup+(-1+h)*ssup*ssr+(-1+h*ssr)*ssrsup)*ssry+(-1+
          h*ssr)*ssry^2),
          1+((ssup^2-2*ssup*ssr+2*h*ssup*ssr-
          h*ssup^2*ssr+4*dsup^2*(-1+h*ssr)*(-1+ssrsup)^2-2*ssup*ssrsup+
          2*h*ssr*ssrsup+2*ssup*ssr*ssrsup-2*h*ssup*ssr*ssrsup+ssrsup^2-
          h*ssr*ssrsup^2-
          2*d*(ssup+h*ssup*ssr+h*ssr*(-2+ssrsup)-ssrsup)*(-1+
          ssry)+(ssup+h*ssup*ssr+h*ssr*(-2+ssrsup)-ssrsup)*ssry-
          2*dsup*(-1+ssrsup)*(-ssr*(2*ssup+h*(-2+ssry))+
          2*d*(-1+h*ssr)*(-1+ssry)+ssry))/(-ssup^2+h*ssup^2*ssr-
          4*dsup^2*(-1+h*ssr)*(-1+ssrsup)^2+2*ssup*ssrsup-2*ssup*ssr*ssrsup+
          2*h*ssup*ssr*ssrsup-ssrsup^2+h*ssr*ssrsup^2+
          4*dsup*(-1+ssrsup)*(-ssup*ssr+2*d*(-1+h*ssr)*(-1+ssry))+
          4*d*ssup*ssr*(-1+ssry)-4*d^2*(-1+h*ssr)*(-1+ssry)^2-
          2*(ssup+(-1+h)*ssup*ssr+(-1+h*ssr)*ssrsup)*ssry+(-1+
          h*ssr)*ssry^2)))
    }
  }
}
#retain those simulations where the suppressor and driver are still segregating
df=subset(df,ysup>=0 & ysup<=1 & xsreqm>=0 & xsreqm<=1)

#Figure A.4 - Equilibrium frequencies of driver and suppressor for a variety of costs in males
ggplot(df,aes(x=xsreqm,y=ysup,col=ssry))+
  geom_point()+ylim(0,1)+xlim(0,0.6)+theme_cowplot()+
  ylab("Equilibrium frequency of Y suppressor")+
  xlab("Equilibrium frequency of sex-ratio X in males")+
  scale_colour_gradient(
     low = "yellow",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour")


#Classify the outcome by Y and X frequency "class"
df$outcome="Neither"
df$outcome[which(df$ysup<0.05)]="Ysup<0.05"
df$outcome[intersect(which(df$xsreqm>0.1),which(df$xsreqm<0.4))]="0.1<Xsr(males)<0.4"
df$outcome[intersect(intersect(which(df$xsreqm>0.1),which(df$xsreqm<0.4)),which(df$ysup<0.05))]="Both"

#Figure A.5 - Conditions for low frequency suppression of drive
ggplot(df,aes(x=ssry,y=ssup,fill=outcome))+geom_tile()+theme_cowplot()+
  xlab("cost of driver in males")+ylab("cost of suppressor")




