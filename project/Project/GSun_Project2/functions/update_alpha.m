function alpha = update_alpha(e,s,C_array,Z_array)
    alpha(1) = (C_array(1)*(-Z_array(1)+Z_array(2)-Z_array(3)));
    alpha(2) = (C_array(2)*(e+Z_array(1)-Z_array(2)+Z_array(3))*(s-Z_array(2)+Z_array(3)));
    alpha(3) = (C_array(3)*(-Z_array(1)+Z_array(2)-Z_array(3)));
%     for alpha_index = 1:length(alpha)
%         alpha(alpha_index) = ...
%             eq(alpha_index,1)*+...
%             eq(alpha_index,2)*+...
%             eq(alpha_index,3)*;
%     end
%     return

% 
%     alpha = C_array(1) * (min([s,(e+Z_array(1))])-Z_array(1));
end