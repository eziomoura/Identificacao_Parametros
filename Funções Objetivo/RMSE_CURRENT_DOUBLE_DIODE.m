%% RMSE (aproximado em uma iter) da corrente - modelo de dois diodos
function [RMSE] = RMSE_CURRENT_DOUBLE_DIODE(x, IVdata, POP_SIZE)
Vmed = IVdata.Vmed(:);
Imed = IVdata.Imed(:);
Vt = IVdata.Vt;
l = length(Imed);
mImed = repmat(Imed.', POP_SIZE,1);
mVmed = repmat(Vmed.', POP_SIZE,1);
Iph = x(:,1);
I01 = x(:,2);
I02 = x(:,3);
n1  = x(:,4);
n2  = x(:,5);
Rs  = x(:,6);
Rp  = x(:,7);
RMSE = sqrt(sum((mImed - Iph + I01.*(exp((mVmed + mImed.*Rs)./(n1*Vt))-1)...
                             + I02.*(exp((mVmed + mImed.*Rs)./(n2*Vt))-1)...
                             + (mVmed + mImed.*Rs)./Rp).^2.')/l);

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