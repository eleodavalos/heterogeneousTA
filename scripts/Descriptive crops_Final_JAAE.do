clear all

global data "C:\Users\zurit\Dropbox\MTEAgricola\Submissions_MTEAgricola\JAAE\Final_Version_JAAE\4github\GithubTorresDavalosMorales2021"
global results "D:\Dropbox\MTEAgricola\Estimaciones\Results"

use "$data\cultivos_y_precios",replace

merge m:1 P_DEPTO  P_MUNIC  UC_UO  ENCUESTA  COD_VEREDA using id_asis
keep if _merge==3
drop _merge

tab P_S11P135,gen(trat)
egen tot_areasem=sum(AREA_SEMBRADA)
gen rend_trat=log(P_S6P57A*Precio*1000/AREA_SEMBRADA) if P_S11P135==1
gen rend_notrat=log(P_S6P57A*Precio*1000/AREA_SEMBRADA) if P_S11P135==2
gen nom_rend=log(P_S6P57A*Precio*1000/AREA_SEMBRADA)
gen real_rend=log(P_S6P57A/AREA_SEMBRADA) if Precio!=.
**Departamentos con precios por cada producto
egen id_dep=tag(P_DEPTO Producto) 
gen precio_dep=Precio if id_dep==1

egen upa=group(UC_UO ENCUESTA COD_VEREDA)

collapse (mean) precio_dep (sum) id_dep AREA_SEMBRADA (mean) tot_areasem (count) upa (sum) trat1 trat2 (mean) real_rend mrend_trat=rend_trat mrend_notrat=rend_notrat (sd) sdrend_trat=rend_trat sdrend_notrat=rend_notrat (count) nrend_trat=rend_trat nrend_notrat=rend_notrat ,by(P_S6P46)
keep if precio_dep!=.
gen per_trat1=trat1/upa
gen per_trat2=trat2/upa
gen per_cropedarea=AREA_SEMBRADA/tot_areasem
gen diff=mrend_trat-mrend_notrat
gen std=((sdrend_trat^2/nrend_trat) + (sdrend_trat^2/nrend_trat))^0.5
gen ttest=abs(diff/std)

keep P_S6P46 precio_dep id_dep per_cropedarea upa per_trat1 per_trat2 diff ttest

export excel using "$results\descriptives_crops.xls",replace firstrow(variables)
