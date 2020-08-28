%% RMSE da corrente - modelo de um diodo
% supondo com T seja desconhecido

function [RMSE] = unknow_T_RMSE_CURRENT_ONE_DIODE(x,IVdata, POP_SIZE)
% const = Ns*k*T/q
Vmed = IVdata.Vmed;
Imed = IVdata.Imed;
const = IVdata.Ns*IVdata.k/IVdata.q;

mImed = repmat(Imed.', POP_SIZE,1);
mVmed = repmat(Vmed.', POP_SIZE,1);
Iph = x(:,1);
I0 = x(:,2);
n = x(:,3);
Rs = x(:,4);
Rp = x(:,5);
Vt = (x(:,6) + 273)*const;
l = length(Imed);
RMSE = sqrt(sum((mImed - Iph + I0.*(exp((mVmed + mImed.*Rs)./(n.*Vt))-1) + (mVmed + mImed.*Rs)./Rp).^2.')/l);
end