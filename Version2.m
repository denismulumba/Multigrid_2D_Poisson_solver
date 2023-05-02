%Convergence plot

C = [84.7777 56.4178 37.3139 24.1225 15.9462 10.4038 6.8261 4.4524 2.9108  1.898...
    1.239 0.80785 0.52703 0.34366 0.22415 0.14617 0.095326 0.062163 0.040539...
    0.026436  0.01724 0.011243 0.0073317 0.0047812 0.0031179 0.0020333 0.001326...
    0.00086469 0.00056389 0.00036773 0.0002398 0.00015638 0.00010198 6.6504e-05 ...
    4.3369e-05 2.8282e-05  1.8444e-05 1.2028e-05 7.8435e-06 5.1149e-06 3.3356e-06...
    2.1752e-06 1.4185e-06 9.2505e-07 6.0325e-07 3.934e-07 2.5654e-07 1.673e-07 ...
    1.091e-07 7.1147e-08];


semilogy(C, 'b-*', LineWidth=2 )
title('Plot to show how the residual norms propergates with the number V-cycles')
xlabel("V-cycles")
ylabel("Residual norm")