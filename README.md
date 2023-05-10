# EnEcon_FSI
Codes to replicate the results in "Time-varying effects of oil price shocks on financial stress: Evidence from India" https://www.sciencedirect.com/science/article/abs/pii/S0140988323002013

This repository contains the following files:
data.csv 
----------------------
Contains the data used to estimate the models
- Date: Period in months
- ROP : month-over-month changes in crude oil prices 
- RXOI:month-over-month changes NYSE Arca Oil Index 
- VIX :Chicago Board Options Exchange's CBOE Volatility Index
- FSI_PCA : Financial Stress Index using Principal Component Analysis
FSIdata.scv
----------------------
contains data used to plot figure 1 and figure 1A in the appendix
FSI_MM : Money market sub-index
FSI_DM : Debt market sub-index
FSI_EM : Equilty market sub-index
FSI_B : Banking  subindex
FSI_FX : Foreign exchange market sub-index
FSI.VEW : Financial Stress Index -variance equal weights (FSI_VEW)
FSI.PCA : Financial Stress Index -Principal comopnents (FSI_PCA)
 
TVPcodes_R2.r
-----------------------
R codes to estimate TVP VAR, to save the IR coefficits as .llo fies and to replicate the figures in the paper. 

TVP_funs_R1.r
-----------------------
Codes for functions used in TVPcodes.r.  

Results_Tables.r
-----------------------
Codes to replicate Table1 and 2 in the paper**
