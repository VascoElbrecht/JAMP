# JAMP
**J**ust **A**nother **M**etabarcoding **P**ipeline - by Vasco Elbrecht - Twitter @luckylionde

JAMP is a modular metabarcoding pipeline, integrating different functions from USEARCH, VSEARCH, CUTADAPT and other programs. The pipeline is run as an R package and automatically generates the needed folders and summary statistics. 

I am still working on the documentation, and the package is in **beta** thus if you run into issues or have questions please contact me per e-mail (luckylion07@googlemail.com).

For a a short tutorial on extracting haplotypes from metabarcoding datasets take a look at the [denoising quick guide](https://github.com/VascoElbrecht/JAMP/wiki/3)-Denoising-quick-guide!).


## Initialling JAMP
Please keep in mind that JAMP needs [Usearch](https://www.drive5.com/usearch/manual/), [Vsearch](https://github.com/torognes/vsearch), and [Cutadapt](cutadapt.readthedocs.io) installed to work properly. Thus Mac or linux based systems are recommended (and windows not officially supported!).

### To install JAMP locally
```# Recommended method
# Installing dependencies needed fro JAMP
install.packages(c("bold", "XML", "seqinr", "devtools", "fastqcr"), dependencies=T)
# Load devtools and install package directly from GitHub
library("devtools")
install_github("VascoElbrecht/JAMP", subdir="JAMP")
```
You can also download the [latest release of JAMP](https://github.com/VascoElbrecht/JAMP/releases), extract and intal within R using `install.packages("JAMP", repos = NULL, type="source")`

### Example of a system wide installation on a ubuntu|debian server:
```bash
wget https://github.com/VascoElbrecht/JAMP/archive/v0.35.tar.gz
tar -xzf v0.35.tar.gz
cd JAMP-0.35
sudo R CMD INSTALL JAMP
```


