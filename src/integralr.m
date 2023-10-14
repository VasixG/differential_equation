function integralr

C = 0.13;

[R_, Rp, Teta_, Tetap] = diffint(2, C*2, 1000, 1000, 0, -10, 10,2000);
time = linspace(-10,10,2000);
%[R_1, Rp1] = diffint(0.03, 1000, 1000, -10, 10,1000);
%[R_2, Rp2] = diffint(0.13, 1000, 1000, -10, 10,1000);
%[R_3, Rp3] = diffint(-0.03, 1000, 1000, -10, 10,1000);
%[R_4, Rp4] = diffint(-0.05, 1000, 1000, -10, 10,1000);
save('T13.mat', 'Teta_');

t = linspace(0,10,1000);
x = zeros(size(t));
y = zeros(size(t));
plot(time, Rp);
hold on;
plot(time, Tetap);
%plot(-F, y);
%plot(F,y);
%plot(t, R_);
%plot(t, Rp4);
%plot(t, Rp1);
%plot(t, Rp2);
%plot(t, Rp3);
xline(0);
% yline(0.785);
yline(0);
%xline(2);
 %legend("R_+(t, (r_*(-1) + r^*(-1))/2, -1)",  "R_+(t, (r_*(-0.05) + r^*(-0.05))/2, -1)", "R_+(t, (r_*(0.03) + r^*(0.03))/2, 0.03)");

for i = 1:length(t)
    x(i) = -Rp(i)*cos(Tetap(i));
    y(i) = -Rp(i)*sin(Tetap(i));
end


hold on;
fimplicit(@(x,y) x.*y-1 - 0.001/2.*exp((x.^2+y.^2)/2), [0 10 0 10]);
plot(x(1),y(1),"*", "MarkerSize",10);
plot(x,y);
hold off;
end

function [R_ , Rp, Teta_, Tetap] = diffint(r0, c, root_t, int_t, t0, time_start, time_end, time_t)
    %r0 = 1;
    f =@(x)1./sqrt(x.^2 - x.^(-2).*(c.*exp(x.^2/2)+2).^2);
    roots1 = @(r) (r^2-2)-c*exp(r^2/2);
    roots2 = @(r) (-r^2-2)-c*exp(r^2/2);

    h = @(r)r.^(-2) * (c*exp(r.^2/2) + 2);

    psi_p = @(r) -1/2 * asin(h(r)) + pi/2;
    psi_m = @(r) 1/2 *asin(h(r));

    list_ofroots1 = [];
    list_ofroots2re = [];
    list_ofroots1re = [];
    list_ofroots2 = [];
    x = linspace(0,20, root_t);
    if c > 0
        for i = x
            val = fzero(roots1, i);
            if ~ ismember( round(val,4) , list_ofroots1) && ~isnan(val) && val >0
                list_ofroots1 = [list_ofroots1, round(val,4)];   
                list_ofroots1re = [list_ofroots1re, val];
            end
    
        end
    else
        for i = x
            val = fzero(roots1, i);
            if ~ ismember( round(val,4) , list_ofroots1) && ~isnan(val) && val > 0
                list_ofroots1 = [list_ofroots1, round(val,4)];   
                list_ofroots1re = [list_ofroots1re, val];
            end
    
            val2 = fzero(roots2, i);
            if ~ ismember( round(val2,4) , list_ofroots2) && ~isnan(val2) && val2 >0
                list_ofroots2 = [list_ofroots2, round(val2,4)];   
                list_ofroots2re = [list_ofroots2re, val2];
            end
    
        end
    end
    
    if c < 0
        if length(list_ofroots1re) == 1 && length(list_ofroots2re) == 1

            a = list_ofroots1re;
            b = list_ofroots2re;

            r0 = a;

            y = linspace(a,b, int_t);

            F = zeros(size(y));

            omega = integral(f,a,b);

            t0_ = integral(f,a,r0);
            t0p = integral(f,r0,b);

            for i = 1:length(y)
                F(i) = real(integral(f,r0,y(i)));
            end

            F_ = -F;

            t = linspace(time_start,time_end,time_t);

            R_ = zeros(size(t));
            Rp = zeros(size(t));

            Teta_ = zeros(size(t));
            Tetap = zeros(size(t));

            for j = 1:length(t)
                k2 = floor((t(j)-t0_-t0)/omega);
                if mod(k2,2) == 0
                    R_(j) = interp1(F, y, t(j)-2*t0_-k2*omega-t0, "linear");
                    if mod(k2, 4) == 0
                        Teta_(j) = psi_m(R_(j));
                    end

                    if mod(k2, 4) == 2
                        Teta_(j) = psi_m(R_(j))+pi;
                    end
                else
                    R_(j) = interp1(F_, y, t(j)-(k2+1)*omega-t0, "linear");
                    if mod(k2, 4) == 1
                        Teta_(j) = psi_p(R_(j))+pi;
                    end
                    if mod(k2, 4) == 3
                        Teta_(j) = psi_p(R_(j));
                    end
                end
    
                k2p = floor((t(j)+t0_-t0)/omega);
    
                if mod(k2p,2) == 1
                    Rp(j) = interp1(F_, y, t(j)-2*t0p-(k2p-1)*omega - t0, "linear");
                    if mod(k2p, 4) == 1
                        Tetap(j) = psi_p(Rp(j))+pi;
                    end

                    if mod(k2p, 4) == 3
                        Tetap(j) = psi_p(Rp(j));
                    end
                else
                    Rp(j) = interp1(F, y, t(j)-(k2p)*omega -t0, "linear");
                    if mod(k2p, 4) == 0
                        Tetap(j) = psi_m(Rp(j));
                    end

                    if mod(k2p, 4) == 2
                        Tetap(j) = psi_m(Rp(j))+pi;
                    end
                end
    
            end
    
        end
    else
        if length(list_ofroots1re) == 2
            a = list_ofroots1re(1);
            b = list_ofroots1re(2);
            
            r0 = a;
            y = linspace(a,b, int_t);
            F = zeros(size(y));
            omega = real(integral(f,a,b));
            t0_ = real(integral(f,a,r0));
            t0p = real(integral(f,r0,b));
    
            for i = 1:length(y)
                F(i) = real(integral(f,r0,y(i)));
            end
    
            F_ = -F;
            t = linspace(time_start,time_end,time_t);
            R_ = zeros(size(t));
            Rp = zeros(size(t));
            Teta_ = zeros(size(t));
            Tetap = zeros(size(t));
            for j = 1:length(t)
                k2 = floor((t(j)-t0_-t0)/omega);
                if mod(k2,2) == 0
                    R_(j) = interp1(F, y, t(j)-2*t0_-k2*omega-t0, "linear");
                    Teta_(j) = psi_m(R_(j));
                else
                    R_(j) = interp1(F_, y, t(j)-(k2+1)*omega-t0, "linear");
                    Teta_(j) = psi_p(R_(j));
                end
    
                k2p = floor((t(j)+t0_-t0)/omega);
    
                if mod(k2p,2) == 1
                    Rp(j) = interp1(F_, y, t(j)-2*t0p-(k2p-1)*omega-t0, "linear");
                    Tetap(j) = psi_p(Rp(j));
                else
                    Rp(j) = interp1(F, y, t(j)-(k2p)*omega-t0, "linear");
                    Tetap(j) = psi_m(Rp(j));
                end
    
            end
        end
    end
    text = sprintf("Для константы c = %f: r_* = %f, r^*  = %f", c/2, a, b);
    display(text);
end
