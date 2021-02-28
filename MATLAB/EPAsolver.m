function [r_s, t_s, v_min_s, v_max_s,isInequalityMet] = EPAsolver(n, alpha, F, display)
    % the solver for equal priority auction
    % n       : the number of buyers;          
    %           should > 1
    % alpha   : probability that the buyer is uninformed;
    %           should be in the range (0,1)
    % F       : the CDF of values;
    %           should be F(0) -> 0; F(1) -> 1
    % display : display solutions if true;
    %           default is true
    
    syms r v_min v_max t 
    syms s v p;
    F(v) = F;
    f = diff(F);
    pi(p) = (1-F(p))*p;
    phi(v) = v - (1-F(v))/f;
    chi(v_min, v_max) = (((1-alpha)*F(v_max)+alpha)^n-((1-alpha)*F(v_min))^n) / (n*(1-alpha)*(F(v_max)-F(v_min))+alpha);

    eq1 = (alpha * (pi(t)-phi(v_max)) == (1-alpha)*((v_min-t)*(phi(v_max)-phi(v_min))*f(v_min)+(F(v_max)-F(v_min))*phi(v_max)-(pi(v_min)-pi(v_max))));
    eq2 = (-alpha*diff(pi(t),t) == (1-alpha)*(phi(v_max)-phi(v_min))*f(v_min));
    eq3 = (-phi(r)*f(r) == (phi(v_max)-phi(v_min))*f(v_min));
    eq4 = ((1-alpha)^(n-1) * int(F(s)^(n-1),s,r,v_min) == chi(v_min,v_max)*(v_min-t));

    [r_s, t_s, v_min_s, v_max_s] = vpasolve([eq1 eq2 eq3 eq4],[r t v_min v_max],[[0,1]; [0,1]; [0,1]; [0,1]]);


    
    % cast output to doubles
    r_s = double(r_s);
    t_s = double(t_s);
    v_min_s = double(v_min_s);
    v_max_s = double(v_max_s);

    % check the inequality conditon
    if (v_max_s >= v_min_s && v_min_s >= t_s && t_s >= r_s)
        isInequalityMet = true;
    else
        isInequalityMet = false;
    end
    
    % whether to display output, default is not to display
    if (nargin == 3)
        display = true;
    end
    if display
        fprintf("r     : %.5f\n", r_s);
        fprintf("t     : %.5f\n", t_s);
        fprintf("v_min : %.5f\n", v_min_s);
        fprintf("v_max : %.5f\n", v_max_s);
        if isInequalityMet
            fprintf("The inequality conditions are met.\n");
        else
            fprintf("WARNING: The inequality conditions are NOT met.\n");
        end
    end
end



