
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Multisource Graph Synthesis with EHR data

[![CRAN](https://www.r-pkg.org/badges/version/PheCAP)](https://CRAN.R-project.org/package=PheCAP)
—– This needs to be changed once we have the links

## Overview

# MUGS

We develop MUlti-source Graph Synthesis (MUGS), an algorithm designed to
create embeddings for pediatric EHR codes by leveraging graphical
information from three distinct sources: (1) pediatric EHR data, (2) EHR
data from the general patient population, and (3) existing hierarchical
medical ontology knowledge shared across different patient populations.

Utilizing existing hierarchical medical ontology as prior general
knowledge, MUGS facilitates efficient transfer learning by grouping
similar codes, thereby enhancing the transferability of knowledge from
general to pediatric systems. To address the heterogeneity within code
groups and between sites, we propose to decompose a code embedding into
three components: the group effect, defined based on the hierarchical
medical ontology; the site-nonspecific code effect, capturing
characteristics of a code that differ from its group effect and are
shared between health systems; and the code-site effect, identifying
site-specific traits of a code. Importantly, this decomposition, coupled
with penalty functions applied to the code and code-site effects,
provides adaptability to varying degrees of heterogeneity within code
groups and between sites and protects against the adverse effects of
negative knowledge transfer via hyperparameter tuning.

<figure>
<img src="man/figures/MUGSFlowchart.png" alt="Flowchart" />
<figcaption aria-hidden="true">Flowchart</figcaption>
</figure>

First, we obtain the two sets of initial embeddings at each site by
using ‘get_embed’ function, then align them by solving the orthogonal
procrustes problem. Second, we utilize the existing hierarchical medical
ontology of PheCodes, LOINC codes, and RxNorms to group similar codes
(<https://shiny.parse-health.org/hierarchies/>), and train initial
embeddings for group effects, code effects, and code-site effects by
pooling the two sets of aligned embeddings. Third, we commence our core
algorithm: updating group effects, code effects, and code-site effects
in an alternating and iterative fashion. ‘GroupEff_par’ and
‘CodeSiteEff_l2_par’ are used to update group effects and code-site
effects, respectively, utilizing parallel computations across multiple
cores or machines to enhance speed. ‘CodeEff_Matrix’ is used to update
code effects via matrix computations.

For hyperparameter tuning, we leverage code-code pairs curated from
related literature. This helps us select the optimal tuning parameters
associated with the penalties on code effects and code-site effects
without the need for data splitting. The performance of different sets
of embeddings with different tuning parameters is evaluated using
‘evaluation.sim’. It is designed to assess the accuracy of the
embeddings in identifying established related pairs versus random pairs
across a wide range of settings.

Although we cannot share the real data from MGB and BCH used to generate
our MUGS embeddings, we have developed a Shiny App
(<https://shiny.parse-health.org/multi-view-net/>) to support downstream
tasks such as pediatric feature engineering, the construction of
pediatric knowledge graphs, and more.

## Installation

Install stable version from CRAN:

``` r
install.packages("MUGS")
```

Install development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("Xuanmengcai/MUGS")
```

## Getting Started

Load in the simulated data (…) and try the R codes from the vignette
\[MUGS.Rmd\] (…).

To use real EHR data, first convert the data into the same format as the
simulated data. Detailed guidelines for data formatting are included in
MUGS.Rmd. Once formatted, call the main function ‘MUGS’ with your data
as input.

## Citations

Li, M., Li, X., Pan, K., Geva, A., Yang, D., Sweet, S. M., Bonzel,
C.-L., Panickan, V. A., Xiong, X., Mandl, K., & Cai, T. (2024).
Multisource representation learning for pediatric knowledge extraction
from electronic health records. npj Digital Medicine.
<https://doi.org/10.1038/s41746-024-01320-4>
