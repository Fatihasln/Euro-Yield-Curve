Usage
r# Step 1: Data exploration & stationarity tests
source("src/01_data_exploration.R")

# Step 2: Extract Nelson-Siegel factors (optional analysis)
source("src/02_dns_extraction.R")

# Step 3: Run forecasting models
source("src/03_forecasting.R")
Note: Place ECB data files in data/raw/ before running.

Results
Performance (Out-of-Sample, 20% Test Set)
MaturityBest ModelRMSFEMAFEImprovement vs RW3MLASSO-VAR0.03300.02290.3%1YLASSO-VAR0.04110.02691.2%5YLASSO-VAR0.05770.04170.0%10YRW0.05550.0414-30YRW0.05350.0404-
Key Findings:

All models perform similarly, confirming Euro sovereign bond market efficiency
LASSO regularization provides marginal improvements (0-1.2%) for short-term maturities
Random Walk remains unbeaten for long-term forecasts (10Y, 30Y)
Results support the "random walk hypothesis" in efficient financial markets


Methodology
Models Implemented

Random Walk (RW) - Benchmark

   ΔY_{t+1} = μ + ε_{t+1}

ARIMA(p,d,q) - Benchmark

Automatic order selection via AIC/BIC
Varying orders by maturity (3M: (5,0,4), others: (0,0,1))


Vector Autoregression (VAR)

   Y_t = Φ₁Y_{t-1} + ... + ΦₚY_{t-p} + ε_t

Captures cross-maturity dynamics
BIC-based lag selection (max lag: 5)


LASSO-VAR

Penalized VAR with L1 regularization
Cross-validated penalty parameter
Prevents overfitting in high-dimensional settings



Estimation Strategy

Rolling fixed window: 80% train, 20% test
1-step ahead forecasts: Out-of-sample evaluation
Model selection: BIC for lag order, CV for LASSO penalty


Key Findings

Stationarity: Yield levels non-stationary (ADF p > 0.05), first differences stationary (p < 0.01)
Market Efficiency: Minimal predictability beyond random walk benchmark
LASSO Benefits: Marginal improvements for short-term maturities through regularization
Cross-Sectional Dynamics: VAR models capture yield curve co-movements but do not significantly improve forecasts


Project Structure
Euro-Yield-Forecasting/
├── data/
│   └── raw/              # ECB CSV files (not included)
├── src/
│   ├── 01_data_exploration.R
│   ├── 02_dns_extraction.R
│   └── 03_forecasting.R
├── results/
│   ├── figures/          # Generated plots
│   └── tables/           # Performance metrics
└── docs/                 # Documentation

Data Source
Provider: European Central Bank (ECB)
Dataset: Euro Area AAA-rated sovereign yield curves
Frequency: Daily (aggregated to quarterly)
Period: 2004-2024
Access: ECB Statistical Data Warehouse

References

Diebold, F. X., & Li, C. (2006). "Forecasting the term structure of government bond yields." Journal of Econometrics, 130(2), 337-364.
Nelson, C. R., & Siegel, A. F. (1987). "Parsimonious modeling of yield curves." Journal of Business, 60(4), 473-489.


Author
Fatih Çiftaslan
Master's Student in Data Science and Economics
University of Milan
📧 mehmetfatih.ciftaslan@studenti.unimi.it
🔗 LinkedIn

License
MIT License - See LICENSE file for details.
This project is for academic purposes. Data belongs to the European Central Bank.

Acknowledgments

European Central Bank for providing high-quality yield curve data
University of Milan, Department of Economics

