# Multiomics Integration for Colorectal Cancer Survival Prediction

[![R](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## 📋 Overview

Machine learning pipeline integrating **clinical, transcriptomic, proteomic, and imaging data** from TCGA to predict colorectal cancer survival outcomes. This work demonstrates how multiomics integration significantly enhances survival prediction accuracy, achieving a **26% improvement** over clinical variables alone.

### Key Achievements
- **C-index: 0.78** (ensemble model)
- **External validation: 0.74** on TCGA-READ dataset
- **Risk stratification**: High-risk patients (2.8y median survival) vs Low-risk (8.2y median survival), p < 0.001
- **Feature reduction**: 530 → 75 features using Elastic Net without performance loss

---

## 🧬 Dataset

**Source**: [UCSC Xena Browser](https://xenabrowser.net/) - TCGA-COAD (Colon Adenocarcinoma)

| Data Type | Samples | Features | Description |
|-----------|---------|----------|-------------|
| **Clinical** | 720 | 10 | Age, stage, grade, MSI status, mutations (KRAS, TP53) |
| **RNA-seq** | 329 | 20,531 → 500 | Gene expression (ANOVA-selected) |
| **RPPA** | 334 | 131 → 20 | Protein expression (Cox-selected: p-AKT, p-ERK, BCL2) |
| **Imaging** | 5,000 | - | CT scans (CNN analysis) |

**Final integrated cohort**: 450 patients with complete multiomics data

---

## 🔬 Methodology

### 1. Data Preprocessing

#### RNA-seq Pipeline
```r
1. Filter low-expression genes (<10% sample presence)
2. Z-score normalization
3. ANOVA F-test for top 500 survival-associated genes
```

#### RPPA Pipeline
```r
1. Univariate Cox regression (131 proteins)
2. Select top 20 proteins (p < 0.05)
3. Focus on pathway targets: AKT, ERK, BCL2, EGFR
```

#### Clinical Features
- 10 key variables: age, stage, grade, tumor size, lymph nodes, metastasis, MSI status, KRAS/TP53 mutations

### 2. Feature Integration

**Elastic Net Regularization** (α = 0.5)
- Input: 500 genes + 20 proteins + 10 clinical = **530 features**
- Output: **75 features** (5-fold CV, lambda.min)

### 3. Model Architecture

#### A. XGBoost Survival Model
```r
Objective: survival:cox
Hyperparameters (Bayesian optimization):
  - eta: 0.01
  - max_depth: 5
  - subsample: 0.8
  - colsample_bytree: 0.8
  - lambda: 1.5, alpha: 0.5

Performance: C-index = 0.82 (5-fold CV)
```

#### B. CNN Imaging Model
- Architecture: ResNet50 + Cox proportional hazards layer
- Input: 5,000 CT scans
- Performance: C-index = 0.65

#### C. Ensemble Model
```r
Final_Score = 0.7 × XGBoost_Score + 0.3 × CNN_Score
Performance: C-index = 0.78
```

---

## 📊 Results

### Model Performance

| Model | Training C-index | Validation C-index | External (TCGA-READ) |
|-------|------------------|-------------------|----------------------|
| Clinical only | 0.62 | 0.58 | 0.56 |
| XGBoost (multiomics) | 0.82 | 0.80 | 0.76 |
| CNN (imaging) | 0.65 | 0.63 | 0.61 |
| **Ensemble** | **0.82** | **0.78** | **0.74** |

### Risk Stratification

| Risk Group | Patients (%) | Median Survival | 5-Year Survival |
|------------|-------------|-----------------|-----------------|
| Low | 33% | 8.2 years | 78% |
| Medium | 33% | 5.1 years | 52% |
| High | 33% | 2.8 years | 24% |

**Log-rank test**: χ² = 42.3, p < 0.001

### Top Predictive Features (SHAP Analysis)

**Genes**: KRAS, TP53, APC, BRAF, PIK3CA  
**Proteins**: p-AKT (Ser473), p-ERK1/2, BCL2, EGFR, p-mTOR  
**Clinical**: Stage, lymph node involvement, MSI status, age

---

## 🛠️ Installation & Usage

### Requirements
```r
# Core packages
install.packages(c("tidyverse", "survival", "glmnet", "caret", "survminer"))

# Machine learning
install.packages(c("xgboost", "ROCR"))

# For imaging (Python required)
pip install tensorflow keras scikit-survival
```

### Run Pipeline
```r
# Load script
source("multiomics_integration_survival.R")

# Execute complete pipeline
results <- main()

# Access components
xgb_model <- results$model
validation_plots <- results$validation$plot
shap_importance <- results$shap
```

---

## 📁 Project Structure

```
colorectal-cancer-multiomics/
│
├── data/
│   ├── TCGA_COAD_clinical.txt
│   ├── TCGA_COAD_rnaseq.txt
│   ├── TCGA_COAD_rppa.txt
│   └── TCGA_COAD_imaging_metadata.csv
│
├── models/
│   ├── final_ensemble_model.rds
│   └── cnn_imaging_predictions.csv
│
├── multiomics_integration_survival.R
├── README.md
└── requirements.txt
```

---

## 🔍 Validation Strategy

1. **Internal validation**: 5-fold cross-validation (C-index: 0.78)
2. **External validation**: TCGA-READ dataset (C-index: 0.74)
3. **Bootstrap resampling**: 1,000 iterations for 95% CI [0.74, 0.82]
4. **Biological validation**: SHAP features align with known CRC pathways

---

## 🎯 Clinical Implications

### Precision Oncology
- **Risk-adapted treatment**: Personalize chemotherapy intensity based on molecular risk
- **Clinical trial stratification**: Enrich high-risk populations for trial enrollment
- **Patient counseling**: Evidence-based prognostic discussions

### Biomarker Discovery
- **Therapeutic targets**: p-AKT, p-ERK for targeted therapy development
- **Prognostic panel**: 75-feature signature for clinical assay development

---

## 📈 Future Directions

- [ ] Multi-center external validation
- [ ] Prospective clinical trial integration
- [ ] Real-time risk calculator web application
- [ ] Integration with treatment response prediction
- [ ] Single-cell RNA-seq for tumor heterogeneity analysis

---

 📚 References

Data Source
- **TCGA-COAD**: The Cancer Genome Atlas - Colon Adenocarcinoma
- **UCSC Xena**: [https://xenabrowser.net/](https://xenabrowser.net/)
Key Publications
1. **TCGA Network** (2012). Comprehensive molecular characterization of human colon and rectal cancer. *Nature*.
2. **Cox Regression**: Cox, D.R. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society*.
3. **XGBoost**: Chen & Guestrin (2016). XGBoost: A Scalable Tree Boosting System. *KDD*.


## 📄 License

MIT License - see [LICENSE](LICENSE) file for details

---

## 🙏 Acknowledgments

- TCGA Research Network for publicly available data
- UCSC Xena Browser team for data curation
- Thesis advisors and collaborators

---

**Note**: This repository contains the analytical pipeline and methodology. Due to data usage agreements, raw TCGA data must be downloaded independently from UCSC Xena Browser.
