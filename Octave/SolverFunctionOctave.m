% standard way of using the function
% EPAsolver(2,0.5,1,true);

% % --- experiment --- %
% % it shows that as alpha increases, v_max would increase, other three would
% % decrease
% for alpha = 0.1:0.1:0.9
%     [r_s, t_s, v_min_s, v_max_s,isInequalityMet] = EPAsolver(3,alpha,0.5);
%     disp([r_s, t_s, v_min_s, v_max_s,isInequalityMet]);
% end
% disp('running');
1;
function [r_s, t_s, v_min_s, v_max_s,isInequalityMet] = EPAsolver(n, alpha, c, display)
    % the solver for equal priority auction
    % n       : should > 1
    % alpha   : should be in the range (0,1)
    % c       : is the exponetial parameter for F(v)i.e. F(v) = v^c,
    % display : display solutions if true; default is false
    
    syms r v_min v_max t 
    syms s v p;
    F(v) = v^(c);
    f = diff(F);
    pi(p) = (1-F(p))*p;
    phi(v) = v - (1-F(v))/f;
    chi(v_min, v_max) = (((1-alpha)*F(v_max)+alpha)^n-((1-alpha)*F(v_min))^n) / (n*(1-alpha)*(F(v_max)-F(v_min))+alpha);

    eq1 = (alpha * (pi(t)-phi(v_max)) == (1-alpha)*((v_min-t)*(phi(v_max)-phi(v_min))*f(v_min)+(F(v_max)-F(v_min))*phi(v_max)-(pi(v_min)-pi(v_max))));
    eq2 = (-alpha*diff(pi(t),t) == (1-alpha)*(phi(v_max)-phi(v_min))*f(v_min));
    eq3 = (-phi(r)*f(r) == (phi(v_max)-phi(v_min))*f(v_min));
    eq4 = ((1-alpha)^(n-1) * int(F(s)^(n-1),s,r,v_min) == chi(v_min,v_max)*(v_min-t));

    out = vpasolve([eq1 eq2 eq3 eq4],[r t v_min v_max],[0.5; 0.5; 0.5; 0.5]);
    % cast output to doubles
    r_s = double(out(1));
    t_s = double(out(2));
    v_min_s = double(out(3));
    v_max_s = double(out(4));

    % check the inequality conditon
    if (v_max_s >= v_min_s && v_min_s >= t_s && t_s >= r_s)
        isInequalityMet = true;
    else
        isInequalityMet = false;
    end
    
    % whether to display output, default is not to display
    if (nargin == 3)
        display = false;
    end
    if display
        disp([r_s, t_s, v_min_s, v_max_s,isInequalityMet]);
    end
end


EPAsolver(2,0.5,1,true);

