syms v;

% ----for simple power functions----%
c = 0.4;
F(v) = v^1;
% ----for simple power functions----%


% ----for more complicated functions----%
z = 2000;
F(v) = ((1/2)*(sqrt(z*v^4 + 1) - 1)/v^2)/(1/2*(sqrt(z*1^4 + 1) - 1)/1^2);
% ----for more complicated functions----%

x = 0.01:0.01:(1-0.01);
y = F(x);
plot(x,y);

EPAsolver(3,0.9,F);
