function history = run_simulation(p,e,s,es,c,SA,t_max,iter, run)
    x0 = [p; e; s; es];
    x = x0;
    z = zeros(3,1);
    alpha = zeros(3,1);
    t = 0;
    identity = eye(3);

    history = [t,alpha',z',zeros(1,3),x'];
    while t<t_max(1)
        % Determine alphas
        alpha = update_alpha(x(2),x(3),c,z)';
        % Determine R_ms
        rm = alpha/sum(alpha);
        % Determine reaction quotient
        rxn_q = rxn_number(rm); 
        % Advance DA
        try 
            z = z + identity(:,rxn_q);
        catch 
            fprintf('Rxn terminated early at run %i, iter %i\n',iter,run);
            break;
        end
        % Determine Tau
        tau = exprnd(1/sum(alpha));
        % Advance Tau
        t = t+tau;
        % Advance history % Record data
        x_advance = x0 + SA*z;
        history = [history;t,alpha',z',rm',x_advance'];
    end
end