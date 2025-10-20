% Get stress director from experimentale stress matrix.
%the computation is basically a mathematical transformation of matrix
%diagonalization 
%%%INPUT
%sxx, syy and sxy are the component of the 2D matrix of stress 
%they should be 2d array 
%%%OUTPUT
%Ms : maximum principle stress 
%ms : minimum principle stress 
%xD : in degree, cos of angle of the rotation to get the Ms/ms matrix
%yD : in degree, sin of angle of the rotation to get the Ms/ms matrix
% mu : anisotropic stress 

%last modification : Lucas, 15/10/2025

function [Ms ms xD yD mu sigmaVM] = getstressdirector_inv(sxx, syy, sxy)
    %size of the array 
    nx = size(sxx,2);
    ny = size(sxx,1);
    %initialize results 
    
    Ms = zeros(ny,nx);%hold the maximum principle stress component 
    ms = zeros(ny,nx);%hold the minimum principle stress component 
    xD = zeros(ny,nx); %hold the cos of the main angle direction
    yD = zeros(ny,nx); %hold the sin of the main angle direction
    mu = zeros(ny,nx); %anisotropic stress
    sigmaVM = zeros(ny,nx); %von mise stress

    for i=1:ny
        for j=1:nx

            xx = sxx(i,j);
            yy = syy(i,j);
            xy = sxy(i,j);
            M = [xx, xy;xy, yy]; %local matrix of stress
            [V,D] = eig(M); % diagonalization
            Ms(i,j) = D(1,1); %eigenvalue
            ms(i,j) = D(2,2); %eigenvalue 

       %now we want the highest eigenvalue to be Ms 
       %this is a correction to switch their position if it's not the case

            if Ms(i,j) >= ms(i,j)
                xD(i,j) = V(1,1); %cos theta
                yD(i,j) = V(2,1); % sin theta
            else
                bap = Ms(i,j);
                Ms(i,j) = ms(i,j);
                ms(i,j) = bap;
                xD(i,j) = V(1,2); %-sin theta
                yD(i,j) = V(2,2);%cos theta 
            end

            mu(i,j) = (Ms(i,j) - ms(i,j))/2;
            sigmaVM(i,j) = sqrt(Ms(i,j)^2 + ms(i,j)^2 - Ms(i,j)*ms(i,j));
           
        end
    end
    
end