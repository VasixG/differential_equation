f_2 = @(r)(-1/2 * r^2 - 1)*exp(-r^2/2);
f_1 = @(r)(1/2 * r^2 - 1)*exp(-r^2/2);

fplot(f_1, [0 5]);

hold on;
fplot(f_2, [0 5]);
grid on;
 xline(0);
 yline(0);
 fplot(@(r) exp(-2), [0 6], "--");
 plot(2,exp(-2),'.', 'MarkerSize',10);
 text(2-0.1,exp(-2)+0.1, "(2, exp(-2))", "FontSize",10)
 plot(2,-3*exp(-2),'.', 'MarkerSize',10);
 text(2,-3*exp(-2), "(2, -3exp(-2))", "FontSize",10)
hold off;

lgd = legend("C^*(\rho)", "C_*(\rho)");
lgd.Location = "southeast";
