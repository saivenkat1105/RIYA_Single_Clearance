%Function for signum alng with shaping torque decay functions
function a = usgn(p,t)

%beta =p.beta;
num = max(size(t),[],'all');

%% Hyperbolic Tangent
%if (j == 1)
%     for i = 1 : num
%         for j = 1:size(t,2)
%             if t(i,j)/beta > 0
%                 t(i,j) = tanh(t(i,j)/beta);
%             else
%                 t(i,j) = 0;
%             end
%         end
%     end
%end

%% Sine Wave
%if (j == 2)
%     for i = 1 : num
%         for j = 1:size(t,2)
%             if (t(i,j)/beta > 0 && t(i,j)/beta < pi/(2) )
%                 t(i,j) = sin(t(i,j)/beta);
%             elseif (t(i,j)/beta > 0 && t(i,j)/beta>pi/(2))
%                 t(i,j) = 1;
%             else
%                 t(i,j) = 0;
%             end
%         end
%     end
%end
%% Actual Excitation
% 
%if j ==3
    t(t>=0) = 1;
    t(t<0) = 0;
%end
%%
a = t;
end