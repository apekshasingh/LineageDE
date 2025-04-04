# LineageDE
## Testing Differential Gene Expression Along Lineages (R package)

### Author: Apeksha Singh (apekshasingh@g.ucla.edu)

Gene expression dynamics resulting from pseudotime reconstruction are fit using a negative binomial generalized additive model (NB-GAM).  
Using the NB-GAM, differential expression is determined between different conditions or lineages using a Likelihood Ratio Test.  
LineageDE is implemented to also accept multiple samples of psuedotime input from which labels are permuted to construct a null distribution.

### Installation Instructions:

(if devtools not already installed)  
install.packages("devtools")  
library(devtools)

install_github("apekshasingh/LineageDE")

### Documentation and Examples:

[Link](https://apekshasingh.github.io/LineageDE/)

### Source Code:

[Link](https://github.com/apekshasingh/LineageDE/)

### Publication:
Population dynamics modeling reveals myeloid bias involves both HSC differentiation and progenitor proliferation biases (Singh A, Blood 2025). See supplemental information for more information.
[Link](https://doi.org/10.1182/blood.2024025598)
