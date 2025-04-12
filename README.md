# Predicting and Modeling Recovery Dynamics Post-Acute Kidney Injury: A Multi-Center Study

### Project Overview

Acute Kidney Injury (AKI) significantly impacts hospitalized patients, affecting 10%-25% and leading to substantial morbidity, mortality, and long-term adverse outcomes. This project leverages electronic health records (EHR) data from four healthcare institutions within the Greater Plains Collaborative (GPC) to:

1. Develop predictive models for short-term (within 7 days of initial onset) AKI progression and reversal.
2. Analyze dynamic patient transitions across AKI states, discharge, and death using multi-state survival modeling.

By combining machine learning predictive modeling with multi-state survival analysis, the study provides a robust framework for identifying risk factors and understanding the intricate dynamics of AKI recovery.

### Repository Structure

The repository is organized into two main components:

#### 1. Predictive Modeling (`predictive_modeling/`)

This folder contains code for the predictive modeling component:
- `data_preprocess.ipynb`: Data cleaning and preprocessing.
- `descriptive_analysis.ipynb`: Descriptive statistics and exploratory data analysis.
- `model_train.ipynb`: Model training using CatBoost and regularized logistic regression.
- `model_evaluation.ipynb`: Evaluation of model performance (ROC, PRC, and calibration curves).
- `model_explain.ipynb`: Feature importance and marginal effect analyses based on SHAP value.
- `consort_diagram.R`: R script for generating the CONSORT diagram describing patient inclusion and exclusion criteria.

#### 2. Multi-state Analysis (`multistate_analysis/`)

This folder includes R scripts for multi-state modeling:
- `multistate_analysis_aki1.R`: Multi-state analysis for AKI Stage 1 patients.
- `multistate_analysis_aki23.R`: Multi-state analysis for AKI Stage 2 and 3 patients.

### Requirements

#### Python Environment

Python packages required:
- catboost  
- jupyter  
- matplotlib 
- numpy  
- pandas 
- pickle
- plotly
- scikit-learn  
- scipy  
- seaborn  
- statsmodels

Install via pip:

```bash
pip install catboost jupyter matplotlib numpy pandas pickle plotly scikit-learn scipy seaborn statsmodels
```

#### R Environment

R libraries required:
- arrow
- diagram
- dplyr
- ggplot2
- stringr
- survival
- tidyverse

Install via R:

```R
install.packages(c("arrow", "diagram", "dplyr", "ggplot2", "stringr", "survival", "tidyverse"))
```

### Contact
For questions or further information regarding the code and analyses, please contact:
**Qi Xu**  
Email: [qixu@ufl.edu](mailto:qixu@ufl.edu)

---

We hope this repository will facilitate further research and collaboration in improving clinical outcomes for patients experiencing AKI.
