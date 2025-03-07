%% Use forward/backward/centered differentiation
function [fwd_x1,fwd_x2,bwd_x1,bwd_x2,ctd_x1,ctd_x2] = func_diff(x1,x2,h)

y = (x1^2+x1^2/x2-exp(x2))*log(x1+exp(x2));

%% partial derivative w.r.t x1
% x1ph = x1 plus h
x1ph = x1+h;
y1ph = [(x1ph^2+x1ph^2/x2-exp(x2))*log(x1ph+exp(x2))];
fwd_x1 = (y1ph-y)/h;

% x1mh = x1 minus h
x1mh = x1-h;
y1mh = [(x1mh^2+x1mh^2/x2-exp(x2))*log(x1mh+exp(x2))];
bwd_x1 = -(y1mh-y)/h;

% centered
ctd_x1 = (y1ph-y1mh)/(2*h);

%% partial derivative w.r.t x2
% x2ph = x2 plus h
x2ph = x2+h;
y2ph = [(x1^2+x1^2/x2ph-exp(x2ph))*log(x1+exp(x2ph))];
fwd_x2 = (y2ph-y)/h;

% x2mh = x2 minus h
x2mh = x2-h;
y2mh = [(x1^2+x1^2/x2mh-exp(x2mh))*log(x1+exp(x2mh))];
bwd_x2 = -(y2mh-y)/h;

% centered
ctd_x2 = (y2ph-y2mh)/(2*h);

end
