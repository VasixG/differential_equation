function integralr
    
%[R_, Rp, R, y] = diffint(-1, 1000, 1000, -10, 10,1000);
[R_1, Rp1, R1, y1] = diffint(0.06, 1000, 1000, -4, 4,1000);
[R_2, Rp2, R2, y2] = diffint(0.26, 1000, 1000, -4, 4,1000);
[R_3, Rp3, R3, y3] = diffint(-0.06, 1000, 1000, -4, 4,1000);
[R_4, Rp4, R4, y4] = diffint(2*-0.95, 1000, 1000, -4, 4,1000);
hold on;
t = linspace(-4,4, 1000);
%plot(-F, y);
%plot(F,y);
%plot(R, y);
plot(t,R_4);
plot(t, R_3);
plot(t,R_1);
plot(t,R_2);
save('tf.mat', 't');

save('Rf03.mat', 'R_1');
save('t03.mat', 'y1');

save('Rf13.mat', 'R_2');
save('t13.mat', 'y2');

save('Rfm3.mat', 'R_3');
save('tm3.mat', 'y3');

save('Rfm5.mat', 'R_4');
save('tm5.mat', 'y4');


 xline(0);
 yline(0);
 legend("C = -0.05",  "C = -0.03", "C = 0.03", "C = 0.13");
hold off;
end
function [R_ , Rp, F, y] = diffint(c, root_t, int_t, time_start, time_end, time_t)
    %r0 = 1;
    f =@(x)1./sqrt(x.^2 - x.^(-2).*(c.*exp(x.^2/2)+2).^2);
    roots1 = @(r) (r^2-2)-c*exp(r^2/2);
    roots2 = @(r) (-r^2-2)-c*exp(r^2/2);
    
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
            t0p = real(integral(f,r0,b));
            for i = 1:length(y)
                F(i) = real(integral(f,r0,y(i)));
            end
            F_ = -F;
            t = linspace(time_start,time_end,time_t);
            R_ = zeros(size(t));
            Rp = zeros(size(t));
            for j = 1:length(t)
                k2 = floor((t(j)-t0_)/omega);
                if mod(k2,2) == 0
                    R_(j) = interp1(F, y, t(j)-2*t0_-k2*omega, "linear");
                else
                    R_(j) = interp1(F_, y, t(j)-(k2+1)*omega, "linear");
                end
    
                k2p = floor((t(j)+t0_)/omega);
    
                if mod(k2p,2) == 1
                    Rp(j) = interp1(F_, y, t(j)-2*t0p-(k2p-1)*omega, "linear");
                else
                    Rp(j) = interp1(F, y, t(j)-(k2p)*omega, "linear");
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
            for j = 1:length(t)
                k2 = floor((t(j)-t0_)/omega);
                if mod(k2,2) == 0
                    R_(j) = interp1(F, y, t(j)-2*t0_-k2*omega, "linear");
                else
                    R_(j) = interp1(F_, y, t(j)-(k2+1)*omega, "linear");
                end
    
                k2p = floor((t(j)+t0_)/omega);
    
                if mod(k2p,2) == 1
                    Rp(j) = interp1(F_, y, t(j)-2*t0p-(k2p-1)*omega, "linear");
                else
                    Rp(j) = interp1(F, y, t(j)-(k2p)*omega, "linear");
                end
    
            end
        end
    end
end
