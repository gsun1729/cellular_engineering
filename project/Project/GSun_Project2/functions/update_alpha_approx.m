function alpha = update_alpha_approx(e,s,C_array,Z_array)
    alpha = C_array(1) * (min([s,(e+Z_array(1))])-Z_array(1));
end