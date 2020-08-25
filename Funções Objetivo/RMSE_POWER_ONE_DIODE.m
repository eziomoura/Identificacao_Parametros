%% RMSE da potência, modelo de um diodo

% a fazer

function [RMSE] = RMSE_POWER_ONE_DIODE(x, Vmed, Imed, Vt, POP_SIZE)
% Imed and Vmed must be a collum vector
l = length(Imed);
Pmed = Vmed.*Imed;
% RMSE = sqrt(sum((mPmed - Iph + I0.*(exp((mVmed + mImed.*Rs)./(n*Vt))-1) + (mVmed + mImed.*Rs)./Rp).^2.')/l);

RMSE = zeros(1,POP_SIZE); % pre-alocacao de memoria
l = length(Imed);
for i = 1:POP_SIZE
    Iph = x(i,1);
    I0 = x(i,2);
    n = x(i,3);
    Rs = x(i,4);
    Rp = x(i,5);
    RMSE(i) = sqrt(sum((Pmed - Vmed.*Iph + Vmed.*I0.*(exp((Vmed+Imed*Rs)/(n*Vt))-1) + Vmed.*(Vmed+Imed*Rs)/Rp).^2)/l);
end
end