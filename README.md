# JAMP
**J**ust **A**nother **M**etabarcoding **P**ipeline - by Vasco Elbrecht - Twitter @luckylionde

JAMP is a modular metabarcoding pipeline, integrating different functions from USEARCH, VSEARCH, CUTADAPT and other programs. The pipeline is run as an R package and automatically generates the needed folders and summary statistics. 

I am still working on the documentation, and the package is in **beta** thus if you run into issues or have questions please contact me per e-mail (luckylion07@googlemail.com).

For a a short tutorial on extracting haplotypes from metabarcoding datasets take a look at the [denoising quick guide](https://github.com/VascoElbrecht/JAMP/wiki/3)-Denoising-quick-guide!).


## To install JAMP locally
Make sure the packages needed for JAMP are installed in R by running `install.packages(c("bold", "XML", "seqinr"), dependencies=T)`. Downlaod the [latest release of JAMP](https://github.com/VascoElbrecht/JAMP/releases), extract and intall within R using `install.packages("JAMP", repos = NULL, type="source")`

## Example of a system wide installation on a ubuntu|debian server:
```bash
wget https://github.com/VascoElbrecht/JAMP/archive/v0.35.tar.gz
tar -xzf v0.35.tar.gz
cd JAMP-0.35
sudo R CMD INSTALL JAMP
```


