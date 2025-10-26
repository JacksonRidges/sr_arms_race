#Load required libraries
library(metR)
library(ggplot2)

#Set up data frame for output (see readme for parameters)
results=data.frame(ssry=numeric(),ssup=numeric(),t=numeric(),SRf=numeric(),SRm=numeric(),Ysup=numeric())
ycrit=0.05 #critical Y(sup) frequency for "invasion"

for(ssry in seq(0,0.4,0.01)) { #loop over ssrys
for(ssup in seq(0,1,0.01)) { #loop over ssups
dY=1/2; dsup=0; ssr=0.59; h=0.135; #define drive strength for naive Y (dY), suppressor Y (dsup); cost in females and dominance of that cost - based on literature values
ssrsup=1-(1-ssry)*(1-ssup) #cost of carrying both driver and suppressor which is just the multiplicative fitness of the two

xf=0.67; xm=0.67; y=0.999 #Starting frequency of drive ~0.33 and Y(sup) = 0.001
ystart=y
t=0 #generation zero

while(y>1-ycrit && y<(ystart+1)/2 && t<10000) { #while loop until y reaches ycrit, an intermediate Y or 10000 generations
  #Define average fitness for X in females (wbarxf), X in males (wbarxm) and Y (in males; wbary)
  wbarxf=xf*xm+(1-h*ssr)*(xm*(1-xf)+xf*(1-xm))+(1-ssr)*(1-xf)*(1-xm)
  wbarxm=1/2*xf*y+1/2*(1-ssup)*xf*(1-y)+(1/2+dY)*y*(1-xf)*(1-ssry)+(1/2+dsup)*(1-ssrsup)*(1-y)*(1-xf)
  wbary=1/2*xf*y+1/2*(1-ssup)*xf*(1-y)+(1/2-dY)*y*(1-xf)*(1-ssry)+(1/2-dsup)*(1-ssrsup)*(1-y)*(1-xf)

  #Calculate allele frequencies in the next generation
  xfnext=(xf*xm+1/2*(1-h*ssr)*(xm*(1-xf)+xf*(1-xm)))/wbarxf
  xmnext=(1/2*xf*y+1/2*(1-ssup)*xf*(1-y))/wbarxm
  ynext=(1/2*xf*y+(1/2-dY)*y*(1-xf)*(1-ssry))/wbary
  
  #Redefine allele frequencies and increment generation
  xf=xfnext; xm=xmnext; y=ynext
  t=t+1
}
#append to results data frame, add output as progress indicator
results[length(results$ssup)+1,]=c(ssry,ssup,t,1-xf,1-xm,1-y)
cat(ssry,ssup,t,xf,xm,y,"\n")
}}


#Countour plot of time to Y(crit) in generations - colors are log(10) transformed
ggplot(subset(results,Ysup>ycrit),aes(x=ssup,y=ssry,z=log10(t),fill=log10(t)))+
  geom_tile()+
  theme_bw()+
  xlab("Cost of suppressor in XY(sup) males")+
  ylab("cost of driver in X(SR)Y males")+
  scale_fill_gradientn(colours = terrain.colors(10))+
  geom_contour(col="black",bins=4)+
  metR::geom_text_contour(aes(z = t)) +
  coord_equal()+
  theme(legend.position = c(0.8,0.7))


