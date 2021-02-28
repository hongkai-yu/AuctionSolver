warning('off','all');

syms r t v_min v_max;
syms s v;

alpha = 0.4;
n = 2;
% c = 0.4;

% F(v) = v^c;

z = 2000;
G(v) = (sqrt(z*v^4 + 1) - 1) / v^2;
F(v) = G(v) / G(1);

% F(v) = ((1/2)*(sqrt(z*v^4 + 1) - 1)/v^2)/(1/2*(sqrt(z*1^4 + 1) - 1)/1^2);

f = diff(F);
pi(v) = (1-F(v))*v;
phi(v) = v - (1-F(v))/f;
chi(v_min, v_max) = (((1-alpha)*F(v_max)+alpha)^n-((1-alpha)*F(v_min))^n) / (n*(1-alpha)*(F(v_max)-F(v_min))+alpha);

eq1 = (alpha * (pi(t)-phi(v_max)) == (1-alpha)*((v_min-t)*(phi(v_max)-phi(v_min))*f(v_min)+(F(v_max)-F(v_min))*phi(v_max)-(pi(v_min)-pi(v_max))));
eq2 = (-alpha*diff(pi(t),t) == (1-alpha)*(phi(v_max)-phi(v_min))*f(v_min));
eq3 = (-phi(r)*f(r) == (phi(v_max)-phi(v_min))*f(v_min));
eq4 = ((1-alpha)^(n-1) * int(F(s)^(n-1),s,r,v_min) == chi(v_min,v_max)*(v_min-t));

equations = [eq1 eq2 eq3 eq4];
vars = [r t v_min v_max];

out = vpasolve(equations, vars,[1/2 1/2 1/2 1/2]);

disp(out);
