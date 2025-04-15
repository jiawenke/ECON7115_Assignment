use "/Users/jiawenke/Downloads/CASIF_98_07.dta", clear

gen k = ln(capital+1)
gen l = ln(employee+1)
gen y = ln(output+1)
gen m = ln(inter_input+1)  // intermediate

* Stage-one
gen k2 = k^2
gen l2 = l^2
gen m2 = m^2
gen kl = k*l
gen km = k*m
gen lm = l*m
gen k3 = k^3
gen l3 = l^3
gen m3 = m^3
gen klm = k*l*m

reg y k l kl km lm klm
predict phi_hat

sort firm_code year

* GMM Function
capture program drop acf_moment

program define acf_moment
    syntax varlist [if] [in], at(name)
    tokenize `varlist'
* Residuals
    local y1 `1'  // capital
    local y2 `2'  // labor
   
    
    * Coefficients
    tempname beta_k beta_l beta_m
    scalar `beta_k' = `at'[1,1]
    scalar `beta_l' = `at'[1,2]
    scalar `beta_m' = `at'[1,3]
    
    * temp vars
    tempvar omega omega_lag g_omega xi
    
    * omega(beta_k, beta_l, beta_m)
    qui gen double `omega' = phi_hat - `beta_k'*k - `beta_l'*l - `beta_m'*m
    
    * lagged omega
    qui by firm_code: gen double `omega_lag' = `omega'[_n-1] if year == year[_n-1] + 1
    
    * Markov transition：omega_t = g(omega_{t-1}) + xi_t
    * Approximation
    qui reg `omega' `omega_lag' c.`omega_lag'#c.`omega_lag' if !missing(`omega_lag')
    qui predict double `g_omega' if !missing(`omega_lag')
    qui gen double `xi' = `omega' - `g_omega' if !missing(`omega_lag')
    
    * Moment conditions
    * 1. E[xi_t * k_t] = 0
    * 2. E[xi_t * l_{t-1}] = 0
    replace `y1' = `xi'*k if !missing(`xi')
    replace `y2' = `xi'*l_lag if !missing(`xi')
end

* lagged labor
sort firm_code year
by firm_code: gen l_lag = l[_n-1] if year == year[_n-1] + 1

* GMM
gmm acf_moment, nequations(2) parameters(beta_k beta_l beta_m) instruments(k l_lag) winitial(identity) onestep from(beta_k = 0.7 beta_l = 0.3 beta_m = 0.1)

* package
replace firm_code = subinstr(firm_code, "SZ", "", .)

sort  firm_code  year
gen exit = 0
replace exit=1  if  firm_code[_n]!=firm_code[_n+1]&year!=2007

destring firm_code, gen (firm_code_num)

xtset firm_code_num year
acfest y, free(l) state(k) proxy(m) 
