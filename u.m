%% Function for the stiffness smoothening
function a = u(p,t)

sigma = p.stiffsigma;
num = max(size(t),[],'all');

%% Hyperbolic Tangent
% % 
% for i = 1 : num
%     for j = 1:size(t,2)
%         if t(i,j) > 0
%             t(i,j) = tanh(sigma*t(i,j));
%         else
%             t(i,j) = 0;
%         end
%     end
% end

%% Arc Tangent
% % 
%     for i = 1 : num
%         if t(i) > 0
%             t(i) = 2/pi*atan(sigma*t(i));
%         else
%             t(i) = 0;
%         end
%     end

%% Hyperbolic Cosine

% for i = 1 : num
%     if t(i) > 0
%         t(i) = log(2*cosh(sigma*t(i)))/(sigma);
%     else
%         t(i) = 0;
%     end
% end

%% Step Function
% 
t(t>=0) = 1;
t(t<0) = 0;

end