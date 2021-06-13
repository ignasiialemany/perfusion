function [obj] = setCapillaryDistributions(obj,L_distribution,Z_distribution,varargin)

if strcmp(L_distribution,'Weibull')
    obj.length_func = @(stream) weibull(stream,obj.average_length,obj.sigma_length);
elseif strcmp(L_distribution,'Constant')
    N=varargin{1};
    obj.length_func = @(stream) (obj.inputs.T * obj.inputs.velocity_value)/N;
else
    fprintf("This length distribution has not been implemented yet");
end

if strcmp(Z_distribution,'Anisotropic')
    obj.direction_func = @(stream,mean_direction) vonmises(stream,mean_direction,obj.k_anisotropy);
elseif strcmp(Z_distribution,'Isotropic')
    obj.direction_func = @(x,mean_direction) isotropic(x);
else
    fprintf("This angle distribution has not been implemented yet");
end
            
end

function [direction] = vonmises(stream,mean_direction,k)
%% VON MISES
%Von Mises distribution needs inverse transform
%10**4 points (used in the paper)
[phi_0,theta_0,~] = cart2sph(mean_direction(1),mean_direction(2),mean_direction(3));
theta_0 = pi/2 - theta_0;
R_z = [cos(-phi_0) -sin(-phi_0) 0; sin(-phi_0) cos(-phi_0) 0; 0 0 1];
R_y = [cos(-theta_0) 0 sin(-theta_0) ; 0 1 0; -sin(-theta_0) 0 cos(-theta_0)];
R = R_y*R_z;

theta = linspace(0,pi,10000);
p_theta = (k)/(2*sinh(k)).*exp(k.*cos(theta)).*sin(theta);
cfd_theta = cumtrapz(theta,p_theta);

%Direction in XYZ'
zenith = interp1(cfd_theta,theta,rand(stream));
azimuth = rand(stream)*2*pi;
xp = 1 .* sin(zenith) .* cos(azimuth);
yp = 1 .* sin(zenith) .* sin(azimuth);
zp = 1 .* cos(zenith);

%Directions on cartesian axis XYZ
direction = inv(R) * [xp;yp;zp];
direction = direction' / norm(direction');
end

function [L] = weibull(stream,average,sigma)
%https://journals.ametsoc.org/view/journals/apme/17/3/1520-0450_1978_017_0350_mfewsf_2_0_co_2.xml?tab_body=pdf
    rnd_value = rand(stream);
    x=linspace(0,600,10000);
    k = (sigma/average)^(-1.086); %Shape factor
    c = average/(gamma(1+1/k)); %Scale factor
    cfd = wblcdf(x,c,k);
	L = interp1(cfd,x,rnd_value*0.999);    
end

function [direction] = isotropic(x)
zenith = acos(-1+2*x);
azimuth = 2*pi*x;
xp = 1 .* sin(zenith) .* cos(azimuth);
yp = 1 .* sin(zenith) .* sin(azimuth);
zp = 1 .* cos(zenith);
direction = [xp yp zp];
end