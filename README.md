# GEM-UPF-Analysis

R scripts for analyses of ultra-processed food (UPF) intake, gut microbial genera, and incident Crohn's disease (CD) in the Genetic, Environmental, and Microbial (GEM) Project.

## Overview

This repository contains two complementary analysis scripts:

1. **Microbiome association analysis:** a two-part generalized estimating equation (GEE) approach to assess associations between UPF intake and genus-level microbial detection and abundance.
2. **Time-to-event interaction analysis:** Cox proportional hazards models to assess whether genus-level microbial features modify the association between UPF intake and incident CD.

These are analysis scripts, not a self-contained, executable reproduction package. The participant-level input dataset is **not included** in this repository.

## Repository files

| File | Purpose | Output |
| --- | --- | --- |
| `two_part_adjusted.R` | Two-part, genus-level GEE association analysis with model-component selection and Benjamini–Hochberg adjustment | `NOVA_kathy_kcal_4.0_genus_twopart_adjusted-updated.csv` |
| `Interaction coxph.zero proportion.R` | Genus-level UPF × microbial-feature interaction analyses using Cox models, plus focused Marvinbryantia analyses | `interaction_kathy_kcal4.0_scaled_and_taxa_scaled_or_binarized8-25.csv` |

## Software requirements

The scripts are written in **R** and load the following packages:

- Two-part GEE script: `geepack`, `Matrix`, `aod`, `MESS`, and `doParallel`.
- Cox interaction script: `survival` and `stringr`.

To install the packages used by the scripts:

```r
install.packages(c(
  "geepack", "Matrix", "aod", "MESS", "doParallel",
  "survival", "stringr"
))
```

Exact R/package versions were not recorded in the supplied scripts. The two-part script calls `MESS:::QIC`, an internal package function whose availability or behavior may vary across versions.

## Input data

Both scripts read a locally available CSV file named:

```text
Nova_genus_interaction.csv
```

Key fields referenced in the code include:

| Variable | Role in analysis |
| --- | --- |
| `SubjectID` | Participant identifier; used to derive a family/cluster identifier (`fam`) |
| `kathy_kcal_4.0` | UPF intake exposure measure, using the source dataset's coding and units |
| `Age`, `Sex` | Covariates |
| `Ethnicity`, `Income` | Additional covariates in the two-part GEE analysis |
| `Total.fiber.intake` | Additional covariate in the two-part GEE analysis |
| `total_reads` | Sequencing-depth covariate in the two-part GEE analysis |
| `CD` | Incident CD event indicator; the Cox script codes `CD == 1` as an event |
| `Fellow_up_durations` | Follow-up time variable used in the Cox models (spelling as in the source code) |
| Genus abundance columns | Taxon-specific microbial abundance measurements |

**Important:** Both scripts identify microbial columns using numeric column positions rather than a separately defined taxa list. Input column order must therefore be checked against the original analytical dataset before running the scripts. The two scripts also use slightly different rules to derive `fam` from `SubjectID`; family identifiers should be verified for consistency.

## Analysis 1: Two-part GEE microbiome associations

**Script:** `two_part_adjusted.R`

For each selected microbial genus, the script considers two model components:

- **Detection component:** a binomial-logit GEE for taxon presence (`abundance > 0`) versus absence (`abundance == 0`). The working correlation is set to `independence`.
- **Positive-abundance component:** a Gaussian GEE for the natural logarithm of abundance among samples with nonzero abundance. The working correlation is set to `exchangeable`.

Both components include `kathy_kcal_4.0`, age, sex, sequencing depth, ethnicity, income, and total fiber intake. Family clustering is specified through `fam`, with robust/sandwich standard errors requested in `geeglm`.

The code applies model-fitting safeguards, including a minimum of 10 observations for the positive-abundance component, checks for model-matrix rank deficiency, and restrictions on fitting the binary component for very low or very high zero proportions. Where both components are fitted, the code also calculates a joint Wald test for the UPF coefficients using a block-diagonal covariance matrix.

The script chooses a taxon-level p-value according to the observed proportion of zero abundance values:

- **Zero proportion ≤ 5%:** use the positive-abundance component p-value.
- **Zero proportion ≥ 15%:** use the joint two-part p-value.
- **Zero proportion > 5% and < 15%:** compare the script's QIC values to select the positive-abundance or joint p-value, subject to availability.

Benjamini–Hochberg false discovery rate (BH-FDR) correction is then applied to the selected p-values. The output contains genus labels, zero proportions, component-specific coefficients and p-values, joint p-values, QIC-related fields, and BH-adjusted q-values. Failed or ineligible models can yield missing values.

**Interpretation:** This analysis assesses associations with microbial detection and conditional positive abundance; it does not itself model time to CD diagnosis or establish causality.

## Analysis 2: UPF × microbiome interactions for incident CD

**Script:** `Interaction coxph.zero proportion.R`

The script fits taxon-specific Cox proportional hazards models for incident CD. UPF exposure (`kathy_kcal_4.0`) is standardized using `scale()`. Microbial features are processed according to their zero proportion:

- **More than 15% zeros:** convert abundance to a binary presence/absence indicator.
- **15% zeros or fewer:** standardize the continuous abundance variable using `scale()`.

Each model includes the UPF variable, the microbial feature, their interaction, age, and sex:

```r
Surv(Fellow_up_durations, CD1) ~
  kathy_kcal_4.0_sc * microbial_feature +
  Sex + Age + cluster(fam)
```

The `cluster(fam)` term requests cluster-robust variance estimation for the family grouping. The batch output contains the coefficient name, log hazard ratio (`Estimate`), hazard ratio (`Exp_Coef`), p-value (`Pr_z`), and exponentiated 95% confidence interval for each reported model coefficient. The UPF × microbial-feature term is the interaction of interest; main-effect coefficients should be interpreted conditional on the other model terms.

The script additionally includes exploratory, focused Cox analyses comparing UPF associations among participants classified as *Marvinbryantia* carriers and non-carriers, followed by a corresponding interaction model. These focused analyses depend on the availability and correct construction of `Marvinbryantia_bin`.

**Interpretation:** The Cox analysis evaluates effect modification on the multiplicative hazard-ratio scale. It does not by itself demonstrate biological interaction or causality.

## How to run

1. Install the packages listed above.
2. Make the authorized analytical dataset available locally as `Nova_genus_interaction.csv`. Do not commit restricted participant-level data to a public repository.
3. **Edit the hard-coded `setwd(...)` line** at the start of each script to point to the directory containing the input CSV and where outputs should be written.
4. Confirm that the required variables, taxon-column positions, taxon labels, and derived family IDs match the intended dataset.
5. Run the scripts individually, for example:

```bash
Rscript two_part_adjusted.R
Rscript "Interaction coxph.zero proportion.R"
```

The scripts write their CSV outputs to the working directory. The interaction script also prints focused *Marvinbryantia* model summaries to the R console.

## Repository

https://github.com/MingyueXue-luna/GEM-UPF-Analysis
