syms r v_min v_max t 
syms s v p;

%---------------INPUT WINDOW-------------------%
n = 5;            % should > 1
alpha = 0.2;      % should in (0,1)
z = 2000;
F(v) = ((1/2)*(sqrt(z*v^4 + 1) - 1)/v^2)/(1/2*(sqrt(z*1^4 + 1) - 1)/1^2);
%---------------INPUT WINDOW-------------------%

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

%    [r_s, t_s, v_min_s, v_max_s] = vpasolve([eq1 eq2 eq3 eq4],[r t v_min v_max],[[0,1]; [0,1]; [0,1]; [0,1]]);
out = vpasolve([eq1 eq2 eq3 eq4],[r t v_min v_max],[[0,1]; [0,1]; [0,1]; [0,1]]);
disp(out.r);

% 
% % cast output to doubles
% %    r_s = double(r_s);
% %    t_s = double(t_s);
% %    v_min_s = double(v_min_s);
% %    v_max_s = double(v_max_s);
% r_s = out(1);
% t_s = out(2);
% v_min_s = out(3);
% v_max_s = out(4);
% 
% % check the inequality conditon
% if (v_max_s >= v_min_s && v_min_s >= t_s && t_s >= r_s)
%     isInequalityMet = true;
% else
%     isInequalityMet = false;
% end
% 
% % whether to display output, default is not to display
% % if (nargin == 3)
%     display = true;
% % end
% if display
%     fprintf("r     : %.5f\n", r_s);
%     fprintf("t     : %.5f\n", t_s);
%     fprintf("v_min : %.5f\n", v_min_s);
%     fprintf("v_max : %.5f\n", v_max_s);
%     if isInequalityMet
%         fprintf("The inequality conditions are met.\n");
%     else
%         fprintf("WARNING: The inequality conditions are NOT met.\n");
%     end
% end
