# OpOptimization



## Overview

- **OpOptimization** is an R package built in the **MFIT 5008** course. This package implements the optimization method of the paper Portfolio optimization with covered calls (in the pdf file DiazKwon2019-1).

- The package is designed for optimizing operational processes using advanced mathematical and statistical models. The package provides tools for solving linear programming problems, generating synthetic datasets, and applying optimization algorithms.

---

## Features

In the R package, an Optimization function is defined in the Optimization.R file, which implements the use of Brownian motion to predict asset price changes and uses a linear programming function (lpSolve) to adjust option weights to minimize CVaR.

---

## Parameters

For ease of understanding, the main parameters in the function are listed here. In addition, you can also view a detailed introduction in the introduction document Optimization.Rd file:

-  **S0**:An array recording the initial prices of different assets(eg:S0 <- c(20804.11,5969.34,2406.67))
  
-  **N**:Number of Brownian motion predictions
  
-  **ST**:The end-of-period price of each asset is predicted N times by Brownian motion, with num_w rows and N columns
  
- **alpha**:alpha
 
-  **num_w**:Number of assets
  
-  **num_p**:The number of options on each asset, an array of length num_w
  
-  **C**:The current market price of the option, with num_w rows and sum_p columns
  
-  **K**:{Strike price for each asset, with num_w rows and sum_p columns
  
-  **rf**:{Risk-free rate of return
  
-  **T**:Duration of asset price change research, using days/365
  
-  **Mu**:Each asset Mu in the Brownian motion, an array of length num_w
  
-  **Sigm**:The volatility of each asset in the Brownian motion, an array of length num_w

---

## Data
The data is in **data_df.rda** and **data_k.rda**. data_df corresponds to parameter C, which is the current market price of the option. data_K corresponds to parameter K, which is the option's expiration exercise price

---
## Group members
This project was completed by several members of the Fintech project (in no particular order):

- ZHANG Wei

- YANG Yifan

- SUN Yuanting

- ZHANG Zejing

- LU Yanbo

- CHEN Guilong

