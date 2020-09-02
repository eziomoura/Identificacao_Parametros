%% RMSE (aproximado em uma iter) da corrente - modelo de um diodo
function [RMSE] = RMSE_CURRENT_SINGLE_DIODE(x, IVdata, POP_SIZE)
Vmed = IVdata.Vmed(:);
Imed = IVdata.Imed(:);
Vt = IVdata.Vt;
l = length(Imed);
mImed = repmat(Imed.', POP_SIZE,1);
mVmed = repmat(Vmed.', POP_SIZE,1);
Iph = x(:,1);
I0 = x(:,2);
n = x(:,3);
Rs = x(:,4);
Rp = x(:,5);
RMSE = sqrt(sum((mImed - Iph + I0.*(exp((mVmed + mImed.*Rs)./(n*Vt))-1) + (mVmed + mImed.*Rs)./Rp).^2.')/l);

%     % sqrt foi removido para melhorar performace
% MSE = zeros(1,popSize); % pre-alocacao de memoria
% l = length(Imed);
% for i = 1:popSize
%     Iph = x(i,1);
%     I0 = x(i,2);
%     n = x(i,3);
%     Rs = x(i,4);
%     Rp = x(i,5);
%     MSE(i) = sum((Imed - Iph + I0*(exp((Vmed+Imed*Rs)/(n*Vt))-1) + (Vmed+Imed*Rs)/Rp).^2)/l;
% end
% f = MSE;
end