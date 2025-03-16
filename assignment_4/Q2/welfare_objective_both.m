function neg_welfare = welfare_objective_both(m, tar_1, tar_2, country_index)
    tarr_vec = func_tar(m, tar_1, tar_2);
    [~, ~, ~, ~, welfare] = func_eqm_iter(tarr_vec, m);
    
    if country_index == 1
        neg_welfare = -welfare(1) - welfare(2); % country 1
    elseif country_index == 2
        neg_welfare = -welfare(3) - welfare(4); % country 2
    else
        error('Invalid country index. Use 1 or 2.');
    end
end
