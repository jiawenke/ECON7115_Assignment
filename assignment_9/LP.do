use "/Users/jiawenke/Downloads/CASIF_98_07.dta", clear


** By hand
* first stage
gen k = ln(capital)
gen l = ln(employee)
gen y = ln(output)
gen i = ln(investment)
*LP
* intermediate
gen m = ln(inter_input)

* first stage
gen k2 = k^2
gen m2 = m^2
gen km = k*m
gen k2m = k2*m
gen km2 = k*m2

reg y l k m k2 m2 km
scalar b_l = _b[l]
predict y_hat
gen phi_hat = y_hat - b_l*l

sort firm_code year

capture program drop omega_moment

program define omega_moment
    syntax varlist [if] [in], at(name)
    tokenize `varlist'
    local y1 `1'  // residual #1
    local y2 `2'  // residual #2
    
    * Define Beta
    tempname beta_k beta_m
    scalar `beta_k' = `at'[1,1]
    scalar `beta_m' = `at'[1,2]
    
    * Temp vars
    tempvar omega omega_lag g_omega xi
    
    * omega(beta_k, beta_m)
    qui gen double `omega' = phi_hat - `beta_k'*k - `beta_m'*m
    
    * omega_{n-1}
    qui by firm_code: gen double `omega_lag' = `omega'[_n-1] if year == year[_n-1] + 1
    
    * approximate g(.) = \xi(beta_k, beta_m)
    qui reg `omega' `omega_lag' c.`omega_lag'#c.`omega_lag' if !missing(`omega_lag')
    qui predict double `g_omega' if !missing(`omega_lag')
    qui gen double `xi' = `omega' - `g_omega' if !missing(`omega_lag')
    
    * Moment conditions
    replace `y1' = `xi' if !missing(`xi')  // condition 1: \xi × k=0
    replace `y2' = `xi' if !missing(`xi')  // condition 2: \xi × m=0
end

* GMM
gmm omega_moment, nequations(2) parameters(beta_k beta_m) ///
    instruments(1: k) instruments(2: m) ///
    winitial(identity) twostep ///
    from(beta_k = 0.4 beta_m = 0.6)

* Store
matrix b = e(b)
scalar beta_k_hat = b[1,1]
scalar beta_m_hat = b[1,2]

* Final productivity
gen omega = phi_hat - beta_k_hat*k - beta_m_hat*m
gen tfp = y - beta_k_hat*k - beta_m_hat*m - b_l*l

* package
replace firm_code = subinstr(firm_code, "SZ", "", .)

sort  firm_code  year
gen exit=0
replace exit=1  if  firm_code[_n]!=firm_code[_n+1]&year!=2007

destring firm_code, gen (firm_code_num)

xtset firm_code_num year

prodest y, free(l) proxy(i) state(k) poly(3) method(op)

