##############################
## Extract identifications from an emailed BOLD link
##
## Karl Cottenie
##
## 2018-12-07
##
##############################

library(tidyverse)
library(viridis)
library(rvest)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

# Startup ends here

## Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (dfName), for lists (lsName), for vectors (vcName) (Thanks Jacqueline May)



## copy-paste the link from the BOLD email
bold = read_html("http://boldsystems.org/index.php/IDS_IdentificationRequest/view?jobId=cottenie.identificationRequest_cottenie_0FB8B7FB-BE9F-40D3-AEEC-07809E47214E.H:bold_jobserver-vm:14008747")
bold 

## extract the individual OTU links
boldOtu = bold %>% html_nodes(".ibox-title") %>%  # finds the tokens of the OTUs
  xml_contents() %>%                              # extracts the contents
  as.character() %>%                              # converts it to character
  str_match("IDS_SingleResult(.*)\\\"\\>") %>%    # selects the token information
  .[,2] %>%                                        # selects the second column
  na.omit() %>% # remove the missing values
  # add the IDS to create the full link
  paste("http://boldsystems.org/index.php/IDS_SingleResult", ., sep = "") %>% 
  lapply(function(x){
  read_html(x) %>% html_table() %>% # read the html file and extract the table
    .[[2]] %>% # the second table is the relevant one
    .[-1,] %>% # this is still part of the header row
    as.tibble() %>% 
    mutate(X8 = as.numeric(X8)) %>% # convert the percentages
    rename(X1 = "Phylym", X2 = "Class", X3 = "Order", X4 = "Family", 
           X5 = "Genus", X6 = "Species", X7 = "Subspecies",
           X8 = "Similarity(%)", X9 = "Status")  # rename columns
})

# TODO create a function for this code, that only needs the emailed link as input
# TODO create similar output to Vasco's Bold_web_hack function

