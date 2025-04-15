use "/Users/jiawenke/Downloads/CASIF_98_07.dta", clear

gen k = ln(capital+1)
gen l = ln(employee+1)
gen y = ln(output+1)
gen m = ln(investment+1)  // intermediate
  
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

reg y k l m k2 l2 m2 kl km lm
predict phi_hat

sort firm_code year

* GMM conditions of ACF
capture program drop acf_moment

program define acf_moment
    syntax varlist [if] [in], at(name)
    tokenize `varlist'
  * residuls
    local y1 `1'
    local y2 `2'
    local y3 `3'
    
    * Beta
    tempname beta_k beta_l beta_m
    scalar `beta_k' = `at'[1,1]  // beta_capital
    scalar `beta_l' = `at'[1,2]  // beta_labor
    scalar `beta_m' = `at'[1,3]  // beta_interm
    
    * temp var
    tempvar omega omega_lag g_omega xi
    
    * omega(beta_k, beta_l, beta_m)
    qui gen double `omega' = phi_hat - `beta_k'*k - `beta_l'*l - `beta_m'*m
    
    * lagged omega
    qui by firm_code: gen double `omega_lag' = `omega'[_n-1] if year == year[_n-1] + 1
    
    * approximation & markov process
    qui reg `omega' `omega_lag' c.`omega_lag'#c.`omega_lag' if !missing(`omega_lag')
    qui predict double `g_omega' if !missing(`omega_lag')
    qui gen double `xi' = `omega' - `g_omega' if !missing(`omega_lag')
    
    * Moment conditions
    * 1. E[xi_t * k_t] = 0
    * 2. E[xi_t * l_{t-1}] = 0
    * 3. E[xi_t * m_{t-1}] = 0
    replace `y1' = `xi' if !missing(`xi')
    replace `y2' = `xi' if !missing(`xi')
    replace `y3' = `xi' if !missing(`xi')
end

* lagged labor and interm
sort firm_code year
by firm_code: gen l_lag = l[_n-1] if year == year[_n-1] + 1
by firm_code: gen m_lag = m[_n-1] if year == year[_n-1] + 1

* GMM
gmm acf_moment, nequations(3) parameters(beta_k beta_l beta_m) ///
    instruments(1: k) ///
    instruments(2: l_lag) ///
    instruments(3: m_lag) ///
    winitial(identity) twostep ///
    from(beta_k = 0.3 beta_l = 0.6 beta_m = 0.4)
