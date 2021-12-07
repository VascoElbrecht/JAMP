<a href="#"><img src="images/JAMP_logo.png" alt="JAMP logo" height="150"/><a/>

## JAMP introduction
**J**ust **A**nother **M**etabarcoding **P**ipeline - Twitter: [@VascoElbrecht](https://twitter.com/VascoElbrecht)

JAMP is modular metabarcoding pipeline, integrating different functions from VSEARCH, CUTADAPT and other programs. The pipeline is run as an R package and automatically generates the needed folders and summary statistics. 

**End of 2021 update:** Currently updating the documentation and older functions for a more streamlined experience :smile: .

For a a short tutorial on extracting haplotypes / ESVs from metabarcoding datasets take a look at the [denoising quick guide](https://github.com/VascoElbrecht/JAMP/wiki/3)-Denoising-quick-guide!).


## Initialling JAMP
Please keep in mind that JAMP needs [Vsearch](https://github.com/torognes/vsearch), and [Cutadapt](cutadapt.readthedocs.io) installed to work properly. Thus Mac or linux based systems are recommended (and windows not officially supported, but you can install e.g. a ubuntu shell on your windows system!).

### To install JAMP locally
```# Recommended method
# Installing dependencies needed fro JAMP
install.packages(c("bold", "XML", "seqinr", "devtools", "fastqcr"), dependencies=T)
# Load devtools and install package directly from GitHub
library("devtools")
install_github("VascoElbrecht/PrimerMiner", subdir="PrimerMiner")
install_github("tobiasgf/lulu")
install_github("VascoElbrecht/JAMP", subdir="JAMP")
```
You can also download the [latest release of JAMP](https://github.com/VascoElbrecht/JAMP/releases), extract and intal within R using `install.packages("JAMP", repos = NULL, type="source")`

### Example of a system wide installation on a ubuntu|debian server:
```bash
wget https://github.com/VascoElbrecht/JAMP/archive/v0.53.tar.gz
tar -xzf v0.53.tar.gz
cd JAMP-0.53
sudo R CMD INSTALL JAMP
```

## Licence

**JAMP is for non profit and academic use only.** If you wish to use any aspects of JAMP commercially, please kindly request permission from [Vasco Elbrecht](mailto:luckylion07@googlemail.com) first. Thank you!


