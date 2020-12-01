function alpha = bisect(fun,a,b)
% Finds the split alpha using bisection
TOL = 10;
fl = fun(a);
fr = fun(b);
if(fr*fl>0)
    alpha=10^a;
    return
end
if(abs(fl)<TOL)
    alpha = 10^a;
    return
elseif(abs(fr)<TOL)
    alpha = 10^b;
    return
else
    fc = fun(.5*(a+b));
    if(fr*fc>0)
        b = .5*(a+b);
    else
        a = .5*(a+b);
    end
    alpha = bisect(fun,a,b);
end