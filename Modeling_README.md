# Modeling presented in Appendix 1

### Parameters defined
|parameter/variable|definition|
|----|-------|
|dY|drive strength against susceptible Y chromosome (often 1/2)|
|ysup|drive strength against suppressing Y chromosome (often 0)|
|ssr|cost of drive in homozygous female|
|h|dominance of cost in females (heterozygote fitness is 1-h*ssr)|
|ssup|cost of Y-linked suppressor in males|
|ssrsup|cost of carrying both drive nd Y linked suppressor|
|ssry|cost of carrying the driver in males|
|t|time in generations|
|xf|frequency of wildtype X in females|
|xm|frequency of wildtype X in males|
|y|frequency of non-suppressing Y in males|


### Mathematica Notebook (pseudoobscura_model.nb)
This file has all the initial derivations of equilibria as well as plots for figures A.2 and A.3. It must be run in Mathematica

### Figure A.1 (FigA.1_invastionTime.R)
This is the script that uses the basic model described and determines the time to invasion of Y suppressors. It does not need input, but parameters can be tweaked

### Figures A.4 and A.5 (FigA.4-5.R)
This is the script that uses the equilibrium frequency from the Mathematica notebook to explore how different parameters influence the equilibrium of both X and Y - shown in A.4 and A.5