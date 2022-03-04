
clear all;
clc;


k=[50,10,8,1,0];
[tensor,co1,co2,co3]=computeTensor("FibreDirection_Colinear.mat");
[tensor,e1_50,e2_50,e3_50]=computeTensor("FibreDirection_VonMises_50.mat");
[tensor,e1_10,e2_10,e3_10]=computeTensor("FibreDirection_VonMises_10.mat");
[tensor,e1_8,e2_8,e3_8]=computeTensor("FibreDirection_VonMises_8.mat");
[tensor,e1_1,e2_1,e3_1]=computeTensor("FibreDirection_VonMises_1.mat");
[tensor,iso1,iso2,iso3]=computeTensor("FibreDirection_VonMises_0.0001.mat");

lambda1=[e1_50,e1_10,e1_8,e1_1,iso1];
lambda2=[e2_50,e2_10,e2_8,e2_1,iso2];
lambda3=[e3_50,e3_10,e3_8,e3_1,iso3];

plot(k,lambda1,"*--");
%hold on
%plot(k,lambda2,"*--");
%plot(k,lambda3,"*--");


%createHistogram("STEAM_N105_Anisotropic_100.mat")



function [tensor,alpha,Eigenvalues2,Eigenvalues3] = computeTensor(str)
v=0.5;
b=0.225;
T=1000;
%ms/um^2
load(str);
index_x  = final_pos(:,1) >=0 & final_pos(:,1) <=3000;
index_y  = final_pos(:,2) >=0 & final_pos(:,2) <=3000;
index_z  = final_pos(:,3) >=0 & final_pos(:,3) <=8000;
inside_voxel = index_x & index_y & index_z;
real_phase = - v .* sqrt(b*T) .* nP(inside_voxel,:);
signal_attenuation = abs(mean(exp(-1i*real_phase), 1));
tensor = process_phase(real_phase,b);
[FA,MD,Eval,Evec]=computeParameters(tensor);
alpha=90-rad2deg(atan(abs(Evec(1,1))/abs(Evec(3,1))));
FA_val=FA;
MD_val=MD;
Eigenvalues1=Eval(1);
Eigenvalues2=Eval(2);
Eigenvalues3=Eval(3);
phase_E1 = computePhase(Evec(:,1)',real_phase);
phase_E2 = computePhase(Evec(:,2)',real_phase);
%drawEigenvectors(Evec)
end

function [phase_output] = computePhase(dir,phase)
phase_output = sum(dir .* phase, 2);
end
function [] = drawEigenvectors(Evec)
%plot3([0,0],[0,2],[0,0],"k--");

%plot3([0,0],[0,0],[0,2],"k--");
%plot3([0,2],[0,0],[0,0],"k--");
%set(gca,'visible','off')
quiver3(0,0,0,Evec(1,1),Evec(2,1),Evec(3,1),'linewidth',4,'ShowArrowHead','off');
hold on;
grid on;
quiver3(0,0,0,Evec(1,2),Evec(2,2),Evec(3,2),'linewidth',4,'ShowArrowHead','off');
quiver3(0,0,0,Evec(1,3),Evec(2,3),Evec(3,3),'linewidth',4,'ShowArrowHead','off');
end

function [tensor] = process_phase(phase, bvalue)
% get diffusion tensor

% gradient sampling directions
directions = [1,  1,  0;
              1, -1,  0;
              1,  0,  1;
              1,  0, -1;
              0,  1,  1;
              0,  1, -1];

ndirs = size(directions, 1);
b_matrix = zeros(ndirs, 6);
signal_ratio = zeros(ndirs, 1);
for i = 1:ndirs
    dir = directions(i, :); % [Gx, Gy, Gz]
    % b-matrix (b_xx, b_yy, b_zz, b_xy, b_xz, b_yz)
    b_matrix(i, :) = [dir([1, 2, 3]).^2, 2*dir([1, 1, 2]).*dir([2, 3, 3])] * bvalue;
    % attenuation vector
    phi = sum(dir .* phase, 2); % combine components
    signal_ratio(i) = abs(mean(exp(-1i*phi), 1));
end

% least squares solution
A = b_matrix;
b = log(signal_ratio);
x = -lscov(A, b); % solve

% assign tensor
tensor = zeros(3, 3);
tensor(tril(true(3))) = x([1, 4, 5, 2, 6, 3]); % assign the right spaces
tensor(triu(true(3), +1)) = tensor(tril(true(3), -1)); % mirror over diagonal
end

function [FA,MD,Eval,Evec] = computeParameters(tensor)
[vector, lambda_3x3] = eig(tensor);
lambda = diag(lambda_3x3).';
[~, idx] = sort(lambda, 'descend');
Evec = vector(:, idx); % (E1, E2, E3) eigenvectors
Eval = lambda(:, idx); % (E1, E2, E3) eigenvalues
E1=Eval(1);
MD = trace(lambda_3x3)/3; % mean diffusivity
%theta_angle = rad2deg(atan(Evec(2,3)/Evec(1,3)));

ADC = [lambda(1),lambda(2),lambda(3)];
MD = (lambda(1)+lambda(2)+lambda(3))/3;

FA_second = sqrt(1.5)*norm(ADC-MD)/norm(ADC);

FA = sqrt(1.5)*sqrt(sum((Eval-MD).^2))./sqrt(sum(Eval.^2)); % fractional anisotropy
%                norm(Eval-MD) / norm(Eval)
end


