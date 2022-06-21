function y=ENV_AMP_EST_local(imf)
%%
% this function was created from code contained within the bivariate emd 
% code written by G. Rilling, P. Flandrin, P. Gonçalves and J. M. Lilly,
% "Bivariate Empirical Mode Decomposition", Signal Processing Letters.


[r,c]=size(imf);
NBSYM=2;
% INTERP='cubic';
INTERP='PCHIP';

y=0*imf(1:end-1,:);
for n=1:r
    x=imf(n,:);
    
    t=1:length(x);

    [indmin,indmax,indzer] = locate_extrema(x);


     %boundary conditions for interpolation
		
    [tmin,tmax,xmin,xmax] = boundary_conditions(indmin,indmax,t,x,NBSYM);

% definition of envelopes from interpolation

     envmax = interp1(tmax,xmax,t,INTERP);
     
     y(n,:)=envmax;
     
end
%----------END OF ENV_AMP_EST_local---------------------------------
