function [L] = weibull(rnd_value,average,sigma)
%https://journals.ametsoc.org/view/journals/apme/17/3/1520-0450_1978_017_0350_mfewsf_2_0_co_2.xml?tab_body=pdf
    x=linspace(0,600,10000);
    k = (sigma/average)^(-1.086); %Shape factor
    c = average/(gamma(1+1/k)); %Scale factor
    cfd = wblcdf(x,c,k);
	L = interp1(cfd,x,rnd_value);    
end