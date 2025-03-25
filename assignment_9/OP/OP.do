use " ", clear

* OP

** By hand
* first stage
rename capital k
rename investment i
rename output y
rename employee l
gen k2 = k^2
gen i2 = i^2
gen ki = k*i
gen k2i = k2*i
gen ki2 = k*i2

reg y l k i k2 i2 ki k2i ki2
scalar b_l = _b[l]
predict y_hat
gen phi_hat = y_hat - b_l*l

sort firm_code year

* Define GMM Objective Function
capture program drop omega_moment

program define omega_moment
    syntax varlist [if] [in], at(name)
    local y : word 1 of `varlist'
    
    * Define Beta_k
    tempname beta_k
    scalar `beta_k' = `at'[1,1]
    
    * Some tempvar
    tempvar omega omega_lag g_omega xi
    
    * omega(beta_k)
    qui gen double `omega' = phi_hat - `beta_k'*k
    
    * omega_lag
    qui by firm_code: gen double `omega_lag' = `omega'[_n-1] if year == year[_n-1] + 1
    
    * polynomial approximate g(.) = \xi(beta_k)
    qui reg `omega' `omega_lag' c.`omega_lag'#c.`omega_lag' if !missing(`omega_lag')
    qui predict double `g_omega' if !missing(`omega_lag')
    qui gen double `xi' = `omega' - `g_omega' if !missing(`omega_lag')
    
    * moment condition
    replace `y' = `xi'*k if !missing(`xi')
end

* GMM
gmm omega_moment, nequations(1) parameters(beta_k) instruments(k) twostep from(beta_k = 0.7) 

* Store
matrix b = e(b)
scalar beta_k_hat = b[1,1]

* Final productivity
gen omega = phi_hat - beta_k_hat*k
gen tfp = y - beta_k_hat*k - b_l*l


** By package
replace firm_code = subinstr(firm_code, "SZ", "", .)

sort  firm_code  year
gen exit=0
replace exit=1  if  firm_code[_n]!=firm_code[_n+1]&year!=2007

destring firm_code, gen (firm_code_num)

xtset firm_code_num year

opreg y, exit(exit) state(k) proxy(i) free(l) vce(bootstrap, reps(2))
