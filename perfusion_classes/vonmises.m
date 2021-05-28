function [zenith] = vonmises(value,mean_zenith,k)
%% VON MISES
%Von Mises distribution needs inverse transform
%10**4 points (used in the paper)
x=linspace(-pi,pi,10000);
p=(1/(2*pi*besseli(0,k)))*exp(k*cos(x-mean_zenith));
cfd = cumtrapz(x,p);
zenith = interp1(cfd,x,value);
end