*******Databases*********
global finaldatabase "C:\Users\zurit\Dropbox\MTEAgricola\Submissions_MTEAgricola\JAAE\Final_Version_JAAE\4github_heterogeneousTA"
global results "D:\Dropbox\MTEAgricola\Estimaciones\Results"

cd "$finaldatabase"
*ssc install reghdfe
*ssc install center
*ssc install avar

set more off

use "Base_estimaciones_Final_JAAE",clear

gen log_rend=log_p1-log_areasem1

sum log_rend log_p1 log_areasem1 if trat!=.

***************Globals of control variables and instruments************

*Controles 
global X  tm cons_agro tot_tier P_S11P138
global P P_S15P168 sex_j edu_j edu edad_j P_S15P169 prod P_S15P165
global M area altura dis_capital dis_mayoris dis_bogota dis_otrosmer region_andina region_caribe region_pasifica 

**Instrumentos
global Z des log_cos_mean


********First Stage OLS*************
reghdfe trat $X $P des if log_cos!=. & log_rend!=.,absorb(COD_VEREDA) cl(COD_VEREDA)
gen sample1=e(sample)
test des
estadd scalar F1_fe=r(F)
estimate store fst1_fe

reghdfe trat $X $P log_cos_mean if log_rend!=.,absorb(COD_VEREDA) cl(COD_VEREDA)
gen sample2=e(sample)
test log_cos
estadd scalar F2_fe=r(F)
estimate store fst2_fe

reghdfe trat $X $P $Z if log_rend!=.,absorb(COD_VEREDA) cl(COD_VEREDA)
test $Z
estadd scalar F3_fe=r(F)
estimate store fst3_fe

********Selection Equation (Probit)*************

global V $X $P des log_cos log_cos_mean log_cos_meadian log_cos_weight
foreach x of global V{
bysort COD_VEREDA: egen v_`x'=mean(`x')
}

global V v_tm v_cons_agro v_tot_tier v_P_S11P138  ///
  v_P_S15P168 v_sex_j v_edu_j v_edu v_edad_j v_P_S15P169 v_prod v_P_S15P165


probit trat $X $P $V v_des des if sample1==1, cl(COD_VEREDA)
margins ,dydx($X $P des) post
estimate store se1_fe

probit trat $X $P $V v_log_cos_mean log_cos_mean if sample2==1, cl(COD_VEREDA)
margins ,dydx($X $P log_cos) post
estimate store se2_fe

probit trat $X $P $V v_des v_log_cos_mean $Z if sample2==1, cl(COD_VEREDA)
margins ,dydx($X $P $Z) post
estimate store se3_fe


*****************************Regresion OLS & IV**************************
reghdfe log_rend trat $X $P if log_cos!=.,absorb(COD_VEREDA) cl(COD_VEREDA)
estimate store ols_fe


**********IV regressions***************

**Exposition to conflict
ivreghdfe log_rend $X $P (trat= des) if log_cos!=.,abs(COD_VEREDA) cl(COD_VEREDA) 
estimates store iv1_fe
estadd scalar CV_1 =r(N_clust)

**Only Cost
ivreghdfe log_rend $X $P (trat= log_cos_mean),abs(COD_VEREDA) cl(COD_VEREDA) 
estimates store iv1a_fe
estadd scalar CV_2 =r(N_clust)


**Both instrument
ivreghdfe log_rend $X $P (trat=des log_cos_mean) ,abs(COD_VEREDA) cl(COD_VEREDA)
estimates store iv2_fe
estadd scalar HJ1=e(jp)
estadd scalar CV_3 =r(N_clust)



***********************Export of results at csv format*********************

global labels trat AsisTéncinca tm TenenciaMaquinaria cons_agro ConstruccionesAgro tot_tier ÁreaUPA P_S11P138 TrabajoPermanente fenomenos DesastresNaturalez tenencia TenenciaPrivada raza Predominanciaétnica ///
P_S15P168 Hombresxhogar sex_j Hombresjefesdehogar edu_j Educaciónjefedehogar edu Educaciónhogar edad_j Edadjefedehogar P_S15P169 Edadhogar prod Fuerzadetrabajo P_S15P165 Tamañohogar ///
area Área altura Altura ind_rur Índiceruralidad dis_capital DistanciaCapital dis_mayoris DistanciaMayorista dis_bogota DistanciaBogota dis_otrosmer Distanciaotromercado region_andina RegiónAndina region_caribe RegiónCaribe region_pasifica RegiónPasifica _cons Constante ///
c_trat AsisTéncinca c_log_areasem1 ÁreaSembrada c_tm TenenciaMaquinaria c_cons_agro ConstruccionesAgro c_tot_tier ÁreaUPA c_P_S11P138 TrabajoPermanente c_fenomenos DesastresNaturalez c_tenencia TenenciaPrivada c_raza Predominanciaétnica ///
c_P_S15P168 Hombresxhogar c_sex_j Hombresjefesdehogar c_edu_j Educaciónjefedehogar c_edu Educaciónhogar c_edad_j Edadjefedehogar c_P_S15P169 Edadhogar c_prod Fuerzadetrabajo c_P_S15P165 Tamañohogar ///
c_area Área c_altura Altura c_ind_rur Índiceruralidad c_dis_capital DistanciaCapital c_dis_mayoris DistanciaMayorista c_dis_bogota DistanciaBogota c_dis_otrosmer Distanciaotromercado c_region_andina RegiónAndina c_region_caribe RegiónCaribe c_region_pasifica RegiónPasifica _cons Constante


**OLS and IV
esttab ols_fe iv1_fe iv1a_fe iv2_fe using "$results\olsiv_fe.csv",replace cells(b(star fmt(%9.3f)) se(par)) ///
           stats(N HJ1  CV_1 CV_2 CV_3 r2_a , fmt(%9.0g %9.3f %9.3f %9.3f))      ///
            legend label mtitle("OLS" "Exposure to conflict (ETC)" "Cost Bartik" "ETC+Bartik mean" "ETC+Bartik median" "ETC+Bartik weight" "ETC+Cost") varlabels($labels) star(* 0.1 ** 0.05 *** 0.01)
**First Stage
esttab fst1_fe fst2_fe fst3_fe using "$results\firststage.csv",replace cells(b(star fmt(%9.3f)) se(par)) ///
           stats(N F1_fe F2_fe F3_fe F4_fe F5_fe r2_a, fmt(%9.0g %9.3f %9.3f %9.3f))      ///
            legend label mtitle("1" "2" "3" "4" "5" ) varlabels($labels) star(* 0.1 ** 0.05 *** 0.01)

**Selection Equation (Probit)
esttab se1_fe se2_fe se3_fe using "$results\selectionequation.csv",replace cells(b(star fmt(%9.3f)) se(par)) ///
           stats(N r2_a, fmt(%9.0g %9.3f %9.3f %9.3f))      ///
            legend label mtitle("1" "2" "3" "4" "5" ) varlabels($labels) star(* 0.1 ** 0.05 *** 0.01)
