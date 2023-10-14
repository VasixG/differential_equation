C = -0.03;
cre = [];
omeg = [];

x = linspace(0,10,100);
for c = C
    roots1 = @(r) (r^2-2)-c*exp(r^2/2);
    roots2 = @(r) (-r^2-2)-c*exp(r^2/2);
    F =@(x)1./sqrt(x.^2 - x.^(-2).*(c.*exp(x.^2/2)+2).^2);
    list_ofroots1 = [];
    list_ofroots2re = [];
    list_ofroots1re = [];
    list_ofroots2 = [];
    
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
            omeg = [omeg, 4*integral(F,list_ofroots1re,list_ofroots2re)];
            cre = [cre, c];
        end
    else
        if length(list_ofroots1re) == 2
            omeg = [omeg, 2*integral(F,list_ofroots1re(1),list_ofroots1re(2))];
            display(list_ofroots1re)
            cre = [cre, c];
        end
    end
end
cre = cre/2;

display(omeg);
scatter(cre/2, 2*omeg, 10, 'filled');