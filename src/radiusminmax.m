f1 = @(c,r) r-sqrt(2*c*exp(r^2/2)+2);
f2 = @(c,r) r-sqrt(-(2*c*exp(r^2/2)+2));

x1 = fimplicit(f2, [-1.05 0.18 0 4]).XData;
y1 = fimplicit(f2, [-1.05 0.18 0 4]).YData;

x2 = fimplicit(f1, [-1.05 0.18 0 2]).XData;
y2 = fimplicit(f1, [-1.05 0.18 0 2]).YData;

x3 = fimplicit(f1, [-1.05 0.18 2 4]).XData;
y3 = fimplicit(f1, [-1.05 0.18 2 4]).YData;

display(x3);

save('radx1.mat', 'x1');
save('rady1.mat', 'y1');

save('radx2.mat', 'x2');
save('rady2.mat', 'y2');

save('radx3.mat', 'x3');
save('rady3.mat', 'y3');
