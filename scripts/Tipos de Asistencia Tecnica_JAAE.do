global raw  "C:\Users\zurit\Dropbox\MTEAgricola\Submissions_MTEAgricola\JAAE\Final_Version_JAAE\4github_heterogeneousTA"

clear all
set more off
cd "$data"
**uniendo base upas y produccion
use "$raw\S01_15(Unidad_productora).dta",clear

tab P_S11P135A_SP1 if P_S11P135==1,m

tab P_S11P135A_SP2 if P_S11P135==1,m

tab P_S11P135A_SP3 if P_S11P135==1,m

tab P_S11P135A_SP4 if P_S11P135==1,m

tab P_S11P135A_SP5 if P_S11P135==1,m

tab P_S11P135A_SP6 if P_S11P135==1,m

tab P_S11P135A_SP7 if P_S11P135==1,m

tab P_S11P135A_SP8 if P_S11P135==1,m

tab P_S11P135A_SP9 if P_S11P135==1,m

tab P_S11P135A_SP10 if P_S11P135==1,m
