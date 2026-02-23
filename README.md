# Signaling Entropy Analysis of Osteosarcoma Cell Lines

Network-based signaling entropy rate (SR) computation for in-house osteosarcoma cancer cell lines, following the Teschendorff framework (SCENT). This project integrates transcriptomic data (TPM) with the STRING protein-protein interaction (PPI) network to quantify signaling promiscuity across 68 osteosarcoma samples.

Collaboration between the [Jacob Scott Lab](https://sites.google.com/site/jacobgscott/) at Cleveland Clinic Lerner Research Institute and Case Western Reserve University.

---

## Methodology

### Signaling Entropy Rate

For each sample, we compute the **signaling entropy rate (SR)** — a measure of how uniformly signaling information flows through a protein interaction network given a sample's expression profile.

The pipeline follows Teschendorff et al. (Methods 2014; Nat Commun 2017):

1. **Expression preparation**: TPM values are converted from log-scale to linear scale to ensure all values are strictly positive.

2. **Network integration**: The STRING PPI network (v12.0, high-confidence edges) is mapped from ENSP IDs to HGNC gene symbols via biomaRt, then intersected with the expression matrix. The maximally connected component is extracted.

3. **Stochastic matrix construction via mass action principle**: For each sample with expression vector **x**, edge-level transition probabilities are computed as:

   ```
   p_ij = A_ij * x_j / Σ_k (A_ik * x_k)
   ```

   where A is the binary adjacency matrix. This yields a row-stochastic matrix **P** describing the probability of signaling flow from gene *i* to gene *j*.

4. **Stationary distribution**: Under detailed balance, the invariant measure is:

   ```
   π_i = x_i * Σ_j(A_ij * x_j) / normalization
   ```

5. **Local and global entropy**:
   - Local entropy per gene: `S_i = -Σ_j p_ij * log(p_ij)`
   - Global entropy rate: `SR = Σ_i π_i * S_i`
   - Normalized: `SR_norm = SR / maxSR`, where `maxSR = log(λ_max)` and `λ_max` is the largest eigenvalue of the adjacency matrix.

### Robustness Validation

- **Bootstrap** (100 iterations): Resample genes with replacement from the entropy matrix, compute summary statistics. Tests stability of aggregate entropy to gene composition.

- **Permutation** (100 iterations): Shuffle expression values within each sample column (preserving marginal distribution, breaking gene-label mapping), recompute SR through the same network. Tests whether the network captures real biological signal vs. random.

- **Network perturbation** (100 iterations × 3 rates): Randomly delete and add 5%, 10%, or 20% of edges, recompute SR. Tests sensitivity of entropy to network topology errors.

---

## File Structure

```
├── entropy_matrix.R                        # Steps 1-3: PPI loading, expression loading, ID conversion
├── entropy_corrected.R                     # Step 4: Corrected SR computation (mass action + entropy rate)
├── bootstrap.R                             # Bootstrap robustness analysis
├── permutation_perturbation_corrected.R    # Corrected permutation & perturbation analyses
├── entropy_heatmaps.R                      # Heatmap visualizations
├── permutation_visualizations.R            # Permutation result plots
```

### Key Outputs

| File | Description |
|------|-------------|
| `entropy_matrix_corrected.csv` | Normalized local entropy per gene per sample |
| `SR_per_sample.csv` | Global entropy rate per sample |
| `gene_entropy_summary_corrected.csv` | Per-gene stats (degree, mean entropy, stationary distribution weight) |
| `permutation_corrected_summary.csv` | Permutation z-scores and p-values per sample |
| `perturbation_corrected_rate_*_summary.csv` | Perturbation CV and shift per sample at each rate |
| `bootstrap_cv_per_sample.csv` | Bootstrap coefficient of variation per sample |

---

## Run Order

```r
# 1. Load PPI, expression data, convert IDs
source("entropy_matrix.R")

# 2. Compute corrected entropy rate
source("entropy_corrected.R")

# 3. Bootstrap (uses entropy_matrix_corrected)
source("bootstrap.R")

# 4. Permutation + perturbation (uses adj_m, maxSR, expr_corrected, SR_per_sample)
source("permutation_perturbation_corrected.R")

# 5. Visualizations
source("entropy_heatmaps.R")
source("permutation_visualizations.R")
```

---

## Data Requirements

- **Expression matrix**: TPM matrix (genes × samples) as CSV. Genes as row names, samples as column names.
- **PPI network**: [STRING v12.0](https://string-db.org/) full protein links file (`9606.protein.links.full.v12.0.txt`).

---

## Dependencies

```r
# CRAN
install.packages(c("tidyverse", "igraph", "entropy", "ggplot2", "pheatmap", "RColorBrewer"))

# Bioconductor
BiocManager::install("biomaRt")
```

---

## References

- Teschendorff AE, Sollich P, Kuehn R. *Signalling entropy: A novel network-theoretical framework for systems analysis and interpretation of functional omic data.* Methods. 2014;67(3):282-293. doi:[10.1016/j.ymeth.2014.03.013](https://doi.org/10.1016/j.ymeth.2014.03.013)

- Teschendorff AE, Enver T. *Single-cell entropy for accurate estimation of differentiation potency from a cell's transcriptome.* Nat Commun. 2017;8:15599. doi:[10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599)

- Banerji CRS, Severini S, Caldas C, Teschendorff AE. *Intra-tumour signalling entropy determines clinical outcome in breast and lung cancer.* PLoS Comput Biol. 2015;11(3):e1004115. doi:[10.1371/journal.pcbi.1004115](https://doi.org/10.1371/journal.pcbi.1004115)

- Teschendorff AE, Banerji CRS, Severini S, Kuehn R, Sollich P. *Increased signaling entropy in cancer requires the scale-free property of protein interaction networks.* Sci Rep. 2015;5:9646. doi:[10.1038/srep09646](https://doi.org/10.1038/srep09646)
