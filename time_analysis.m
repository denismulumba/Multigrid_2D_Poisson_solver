

time_Source_term1 = [0.0024259 0.0022816 0.0071199 0.0081227 0.011898 0.031225  0.032831 ...
    0.20029 2.1523 ];
%time_Source_term2 = [0.002478  0.0021696 0.0052827  0.0071663 0.010686 0.025298  0.13084 ...
  %  0.1331  1.1022 ];
T = table(time_Source_term1', 'VariableNames',{'time (S)'});
disp(T)

writetable(T)