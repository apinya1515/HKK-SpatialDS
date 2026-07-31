# Executive Presentation & Manuscript Summary
## Bayesian SECR Distance Sampling Analysis of Ungulate Communities
### Huai Kha Khaeng Wildlife Sanctuary (HKK), Thailand (3,011 km²)

---

## 1. Executive Summary & Key Population Estimates

This study implements a spatially explicit Bayesian SECR distance sampling model incorporating spatial intrinsic CAR (ICAR) random effects and environmental covariates based on **Kumar (2021)**. A total of **5 ungulate species** across **18 top-performing candidate models** ($\Delta \text{wAIC} \le 2.0$) were evaluated.

### Master Ungulate Population & Density Estimates (Study Area = 3,011 km²)

| Species | Species Code | Individual Density ($D_{\text{indiv}} / \text{km}^2$) | Total Abundance ($N_{\text{total}}$) | 95% Bayesian Credible Interval | Cluster Density ($D_{\text{cluster}} / \text{km}^2$) | Average Group Size ($\text{AGS}$) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Muntjac** | `MJK` | **4.69 / km²** | **14,124** | [11,533 - 17,436] | **3.75 / km²** | 1.25 ind/group |
| **Sambar deer** | `SBR` | **3.08 / km²** | **9,286** | [6,523 - 13,288] | **1.71 / km²** | 1.80 ind/group |
| **Wild boar** | `PIG` | **2.12 / km²** | **6,373** | [3,554 - 12,107] | **0.86 / km²** | 2.45 ind/group |
| **Gaur** | `GAR` | **1.10 / km²** | **3,312** | [1,578 - 7,389] | **0.38 / km²** | 2.92 ind/group |
| **Banteng** | `BTG` | **1.07 / km²** | **3,229** | [1,273 - 8,673] | **0.28 / km²** | 3.78 ind/group |

---

## 2. Key Ecological Drivers & Covariate Effects

Across all 5 species, environmental predictors demonstrated clear ecological structuring:

1. **Distance to Streams (`dist_str`)**:
   - **Strong Negative Effect ($\beta < 0$)**: Strongly preferred by **Banteng** ($\beta = -0.582$) and **Sambar deer** ($\beta = -0.412$). Animal densities are significantly higher in close proximity to riparian stream corridors.
2. **NDVI Seasonality / CV (`ndvi_cv`)**:
   - **Positive Effect ($\beta > 0$)**: Selected in top models for **Banteng**, **Gaur**, and **Sambar deer**, indicating preference for edge habitats with high vegetation heterogeneity and forage seasonality.
3. **Elevation (`elev`) & Slope (`slope`)**:
   - **Muntjac & Wild boar**: Positively associated with moderate slope and mid-elevation ridges, avoiding heavily flooded lowlands.
4. **Habitat Cover Types (`BB` Bamboo, `DE` Dry Evergreen, `DD` Deciduous)**:
   - **Gaur**: Strongly associated with **Dry Evergreen (`DE`)** forest types.
   - **Banteng & Sambar**: High inclusion probability for **Bamboo (`BB`)** understory cover.

---

## 3. Methodological Rigor & Diagnostics

- **Gelman-Rubin Diagnostics ($\hat{R}$)**: All primary density and detection parameters converged with $\hat{R} < 1.05$.
- **Effective Sample Size ($\text{ESS}$)**: All model parameters achieved $\text{ESS} > 400$ independent draws.
- **Monte Carlo Standard Error ($\text{MCSE}$)**: $\text{MCSE} < 1\%$ of parameter standard deviation across all 18 models.
- **Multi-collinearity Verification**: Pairwise Pearson correlations among spatial predictors confirmed $|r| < 0.70$.

---

## 4. Key Figures & Presentation Files

1. 📊 **Community Density Comparison**:
   - [Community_Density_Comparison.png](file:///d:/GitHub/HKK-SpatialDS/Results/Density/Community_Density_Comparison.png)
2. 🌡️ **Spatial Covariate Correlation Heatmap**:
   - [Spatial_Covariates_Correlation_Heatmap.png](file:///d:/GitHub/HKK-SpatialDS/Results/Covariates/Spatial_Covariates_Correlation_Heatmap.png)
3. 📉 **3-Panel Density & Group Size Histograms**:
   - [Results/Density/](file:///d:/GitHub/HKK-SpatialDS/Results/Density/)
4. 🗺️ **CAR Spatial Random Effect Maps**:
   - [Results/CAR_map/](file:///d:/GitHub/HKK-SpatialDS/Results/CAR_map/)
5. 📦 **Covariate Beta Boxplots with Quantiles**:
   - [Results/Covariates/](file:///d:/GitHub/HKK-SpatialDS/Results/Covariates/)
6. 📋 **Full Model Selection & Coefficient Master Tables**:
   - [Full_Candidate_Model_Selection_Summary.xlsx](file:///d:/GitHub/HKK-SpatialDS/Results/tables/Full_Candidate_Model_Selection_Summary.xlsx)
   - [Ecological_Covariate_Effects_Synthesis.xlsx](file:///d:/GitHub/HKK-SpatialDS/Results/tables/Ecological_Covariate_Effects_Synthesis.xlsx)
