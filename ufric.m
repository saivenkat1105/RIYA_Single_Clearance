%% Function for friction stiffness
function a = ufric(p,t)

% sigma =p.fricsigma;
sigma = 50;
num = max(size(t),[],'all');

%% Hyperbolic Tangent
 t = tanh(sigma*t);
% for i = 1 : num
%     if t(i) > 0
%         t(i) = tanh(sigma*t(i));
%     else
%         t(i) = 0;
%     end
% end

%% Arc Tangent
% t = 2/pi*atan(sigma*t);
% for i = 1 : num
%     if t(i) > 0
%         t(i) = 2/pi*atan(sigma*t(i));
%     else
%         t(i) = 0;
%     end
% end

%% Hyperbolic Cosine
% t = log(2*cosh(sigma))/(sigma);
% for i = 1 : num
%     if t(i) > 0
%         t(i) = log(2*cosh(sigma*t(i)))/(sigma);
%     else
%         t(i) = 0;
%     end
% end

%% Step Function
% % 
% t(t>0) = 1;
% t(t<0) = -1;
% t(t==0) = 0;
%% Signum HypCos Function
% t = abs(t);
% 
a = t;

end