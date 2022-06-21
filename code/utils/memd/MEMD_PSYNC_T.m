
function [RO_IMFs,T_vec,Pairs]=MEMD_PSYNC_T(Fs,W,W_OV,IMF)

% Edit line 112 % 135 frolm NUM_EXTEMA < 8 to NUM_EXTEMA < 3

% Determines the phase synchrony between IMFs of different channels

% Required functions : 'memd.m' by Rehman and Mandic.

%-----------       OUTPUT     ------------------------------
% 'RO_IMFs' is an M X N X P matrix with each entry corresponding to the 
% the degree of phase synchrony between the same IMF index of two channels,
% where M is the number of all channel pairs, eg. there are 120 channel pairs, given 16 channels (nchoosek(16,2))
% N is the number of IMFs, and P is the number of time windows

% The values are of the range [0 - 1] with 1 indicating max synchrony
% and 0 indicating no sync.
%
%
% 'T_vec' is a 1 X P vector which denotes the time grid for 'RO_IMFs'
%
% 'Pairs' is an M X 2 matrix which shows the order of channel pairs corresponding to M channel pairs of 'RO_IMFs'

% Example:
% Given 10 IMFs of 4 channels obtained using memd.m, the size of 'RO_IMFs' is
% 6 X 10 X P
% 6 is the number of channel pairs determined by nchoosek(4,2)
% 10 is the number of IMFs
%
% 'Pairs' is
%      1     2
%      1     3
%      1     4
%      2     3
%      2     4
%      3     4
% which shows the order of channel pairs correspending to the 'M' of 'RO_IMFs'.


%-----------       INPUT     ------------------------------
% Fs :: Sampling frequency.
% W :: Window length to estimate phase sync.
%      Suggest 200.
% W_OV :: Number of samples for windows overlapping     
%         but if W_OV is a ratio (<1) then W_OV = W_OV * W
% IMF :: IMFs obtained from MEMD


%----------Start of main function---------------------------------

if W_OV<1
    W_shift=1-W_OV;
    W_shift=round(W_shift*W);
else
    W_shift=W-W_OV;
end



len=size(IMF,3);

% Number of channels
num_ch=size(IMF,1);

% All pairs that are available
Pairs=nchoosek(1:num_ch,2);

% Number of pairs
num_pairs=size(Pairs,1);

% Number of IMFs
NUM_IMFs=size(IMF,2);

%%
q_count=0;
for q=1:W_shift:len-W+1
    q_count=q_count+1;
end

for pair_count=1:num_pairs
    
    %% Determing minium index
    imf=zeros(2,size(IMF,2),size(IMF,3));
    imf(1,:,:)=squeeze(IMF(Pairs(pair_count,1),:,:));
    imf(2,:,:)=squeeze(IMF(Pairs(pair_count,2),:,:));
       
%     x1=X(Pairs(pair_count,1),:);
%     x2=X(Pairs(pair_count,2),:);
    
    %- -  Analysis is only performed on IMFs with sufficient numbers of extrema
    % default
    IDX_r=NUM_IMFs;
    % To determine the number of extrema in each IMF
    % for k = 1:NUM_IMFs
    k=1;
    stop=0;
    while ~stop && (k < NUM_IMFs+1)
%           [indmin,indmax] = locate_extrema(real(imf(k,:))); Apit commented out 250414

        % Apit added 250414
        [indmin,indmax] = locate_extrema(squeeze(imf(1,k,:))');

        NUM_EXTEMA(1) = length(indmin) + length(indmax);

        if (NUM_EXTEMA < 3)        
          IDX_r=k-1;
          stop=1;
        end

        k=k+1;
    end


    % default
    IDX_i=NUM_IMFs;
    % To determine the number of extrema in each IMF
    % for k = 1:NUM_IMFs
    k=1;
    stop=0;
    while ~stop && (k < NUM_IMFs+1)
%         [indmin,indmax] = locate_extrema(imag(imf(k,:)));   Apit commented out 250414
        
        % Apit added 250414
        [indmin,indmax] = locate_extrema(squeeze(imf(2,k,:))');
        
        NUM_EXTEMA(1) = length(indmin) + length(indmax);

        if (NUM_EXTEMA < 3)        
          IDX_i=k-1;
          stop=1;
        end

        k=k+1;
    end


    % Analysis index
    IDX=min(IDX_r,IDX_i);


    %% Phase synchrony calculation

    for nn=1:IDX
        
        [ro,Ar,Ai,Fr,Fi,T_vec]=PHASE_SYNC_local(squeeze(imf(1,nn,:))',squeeze(imf(2,nn,:))',Fs,W,W_shift);
        
        RO_IMFs(pair_count,nn,:)=ro; 

    end
        
end

T_vec=(T_vec-1+W_shift)/Fs;
% F_vec=F_vec_norm*Fs;

%----------END OF main function---------------------------------



function [ro,instAmp_r,instAmp_i,instFreq_r,instFreq_i,t_vec]=PHASE_SYNC_local(C1,C2,Fs,W,W_OV)
%%
% returns vetor 'ro' - the degree of phase synchrony between C1(k) and C2(k) (IMFs).
%                 

% IMF1 and IMF2 are i/p vectors
% Fs is the sampling rate
% W is the window length


% W must be a multiple of 2
W=W-rem(W,2);

len=length(C1);
%MEM ALOC.---------------------
 instAmp_r=zeros(1,len);
 instFreq_r=zeros(1,len);
 instAmp_i=zeros(1,len);
 instFreq_i=zeros(1,len);

 % ro=zeros(1,len-(W+1));   Apit commented out 240414
 
 %%%%% Apit added  240414 %%%%%%%%
 q_count=0;
 for q=1:W_OV:len-W+1
     q_count=q_count+1;
 end
 ro=zeros(1,q_count);
 t_vec=zeros(1,q_count);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------

% obtain envelopes for IMFs
y2=abs(ENV_AMP_EST_local(C2));
y1=abs(ENV_AMP_EST_local(C1));

y2=y2+(0.002*mean(abs(y2)));
y1=y1+(0.002*mean(abs(y1)));


% normalise IMFs
zi=C2./y2;
zr=C1./y1;

% obtain instantaneous amplitude, instantaneous frequency and instantaneous
% phase information for normalised IMFs
[instAmp_r_t,instFreq_r_t,PHI_r] = INST_FREQ_local(zr,Fs);
[instAmp_i_t,instFreq_i_t,PHI_i] = INST_FREQ_local(zi,Fs);



instAmp_r=y1;
instFreq_r=(instFreq_r_t);
instAmp_i=y2;
instFreq_i=(instFreq_i_t);



PHI_r=PHI_r+((PHI_r<0).*(-1*PHI_r))+((PHI_r<0).*(1*PHI_r))+pi;
PHI_i=PHI_i+((PHI_i<0).*(-1*PHI_i))+((PHI_i<0).*(1*PHI_i))+pi;

% obtain phase difference between the signals
PHI_diff=zeros(size(PHI_r));
for q=1:length(PHI_r)
    PHI_diff(q)=min(mod(abs(PHI_r(q)-PHI_i(q)),2*pi),...
        mod(abs(PHI_r(q)-PHI_i(q)+(2*pi)),2*pi));
    PHI_diff(q)=min(PHI_diff(q),mod(abs(PHI_r(q)-PHI_i(q)-(2*pi)),2*pi));
end
   

% calculate degree of phase synchrony for a sliding window applied to
% phase difference (PHI_diff)

q_count=0;
  for q=1:W_OV:length(PHI_r)-W+1      % Apit added ':W_OV:' on 240414
    
    L=W;
    M=ceil(exp(.626 + (0.4*log(L))));
    % maximum entropy of phase distribution
    E_max=log(M);
    n_bins=M; %% i.e. n_bins=ceil(exp(E_max));
    edges = linspace(0,1*pi,n_bins);
    h=histc(PHI_diff(q:q+W-1),edges);
    total_num_bin_counts=sum(h);
    % 'probability distribution' for a segment of PHI_diff
    p=h/total_num_bin_counts;
    % entropy of phase distribution
    E=sum(-p.*log(p+(p==0)));
    
    % degree of phase synchrony at location k
    % ro(q)=(E_max-E)/E_max;    Apit commented out 240414
    
    % Apit added 240414
    q_count=q_count+1;
    t_vec(q_count)=q;
    ro(q_count)=(E_max-E)/E_max;
  end    
    
%----------END OF PHASE_SYNC_local---------------------------------

% function [HH_matrix,Fgrid]=HIL_HUANG_SPEC_local(F,A,FR,NFB)
% %%
% % computes Hilbert Huang spectrum for phase synchrony information
% 
% % ouput is HH_matrix an M X N matrix with each entry corresponding to the 
% % the degree of phase synchrony in frequency (M) and time (N)
% 
% % 'Fgrid' defines the frequency components of HH_matrix
% 
% 
% 
% NUM1=((.5-FR)/FR)+1;
% Fgrid=linspace(FR,.5,NUM1)-(FR/2);
% 
% NUM2=((.5)/FR)+1;
% Fgrid_alt=linspace(0, .5 + (FR/2),NUM2);
% 
% 
% 
% if size(F,1)==1
% 
% 
%     for qq=2:NFB+1
%             
%     temp=double((F>=Fgrid_alt(qq-1))&(F<Fgrid_alt(qq)));
%             
%     HH_matrix(qq-1,:)=(A.*temp);
%             
%     end
% end
% 
% 
% if size(F,1)>1
%     for qq=2:length(Fgrid_alt)
%     
%     temp=(double((F>=Fgrid_alt(qq-1))&(F<Fgrid_alt(qq))));
%             
%     HH_matrix(qq-1,:)=sum(A.*temp);
%             
%    end
% end
% 
% 
% 
% %----------END OF HIL_HUANG_SPEC_local---------------------------------


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

function [indmin, indmax, indzer] = locate_extrema(x)
%%
%extracts the indices corresponding to extrema

% this function was created from code contained within the bivariate emd 
% code written by G. Rilling, P. Flandrin, P. Gonçalves and J. M. Lilly,
% "Bivariate Empirical Mode Decomposition", Signal Processing Letters.

t=1:length(x);


m = length(x);

if nargout > 2
	x1=x(1:m-1);
	x2=x(2:m);
	indzer = find(x1.*x2<0);
	
	if any(x == 0)
	  iz = find( x==0 );
	  indz = [];
	  if any(diff(iz)==1)
	    zer = x == 0;
	    dz = diff([0 zer 0]);
	    debz = find(dz == 1);
	    finz = find(dz == -1)-1;
	    indz = round((debz+finz)/2);
	  else
	    indz = iz;
	  end
	  indzer = sort([indzer indz]);
	end
end
  
d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;


% when two or more consecutive points have the same value we consider only 
% one extremum in the middle of the constant area

if any(d==0)
  
  imax = [];
  imin = [];
  
  bad = (d==0);
  dd = diff([0 bad 0]);
  debs = find(dd == 1);
  fins = find(dd == -1);
  if debs(1) == 1
    if length(debs) > 1
      debs = debs(2:end);
      fins = fins(2:end);
    else
      debs = [];
      fins = [];
    end
  end
  if length(debs) > 0
    if fins(end) == m
      if length(debs) > 1
        debs = debs(1:(end-1));
        fins = fins(1:(end-1));

      else
        debs = [];
        fins = [];
      end      
    end
  end
  lc = length(debs);
  if lc > 0
    for k = 1:lc
      if d(debs(k)-1) > 0
        if d(fins(k)) < 0
          imax = [imax round((fins(k)+debs(k))/2)];
        end
      else
        if d(fins(k)) > 0
          imin = [imin round((fins(k)+debs(k))/2)];
        end
      end
    end
  end
  
  if length(imax) > 0
    indmax = sort([indmax imax]);
  end

  if length(imin) > 0
    indmin = sort([indmin imin]);
  end
  
end  
%----------END OF ENV_AMP_EST_local---------------------------------


function [tmin,tmax,xmin,xmax] = boundary_conditions(indmin,indmax,t,x,nbsym)
%%
% computes the boundary conditions for interpolation (mainly mirror symmetry)

% this function was created from code contained within the bivariate emd 
% code written by G. Rilling, P. Flandrin, P. Gonçalves and J. M. Lilly,
% "Bivariate Empirical Mode Decomposition", Signal Processing Letters.

	
	lx = length(x);
	
	if (length(indmin) + length(indmax) < 3)
		error('not enough extrema')
	end

	if indmax(1) < indmin(1)
    	if x(1) > x(indmin(1))
			lmax = fliplr(indmax(2:min(end,nbsym+1)));
			lmin = fliplr(indmin(1:min(end,nbsym)));
			lsym = indmax(1);
		else
			lmax = fliplr(indmax(1:min(end,nbsym)));
			lmin = [fliplr(indmin(1:min(end,nbsym-1))),1];
			lsym = 1;
		end
	else

		if x(1) < x(indmax(1))
			lmax = fliplr(indmax(1:min(end,nbsym)));
			lmin = fliplr(indmin(2:min(end,nbsym+1)));
			lsym = indmin(1);
		else
			lmax = [fliplr(indmax(1:min(end,nbsym-1))),1];
			lmin = fliplr(indmin(1:min(end,nbsym)));
			lsym = 1;
		end
	end
    
	if indmax(end) < indmin(end)
		if x(end) < x(indmax(end))
			rmax = fliplr(indmax(max(end-nbsym+1,1):end));
			rmin = fliplr(indmin(max(end-nbsym,1):end-1));
			rsym = indmin(end);
		else
			rmax = [lx,fliplr(indmax(max(end-nbsym+2,1):end))];
			rmin = fliplr(indmin(max(end-nbsym+1,1):end));
			rsym = lx;
		end
	else
		if x(end) > x(indmin(end))
			rmax = fliplr(indmax(max(end-nbsym,1):end-1));
			rmin = fliplr(indmin(max(end-nbsym+1,1):end));
			rsym = indmax(end);
		else
			rmax = fliplr(indmax(max(end-nbsym+1,1):end));
			rmin = [lx,fliplr(indmin(max(end-nbsym+2,1):end))];
			rsym = lx;
		end
	end
    
	tlmin = 2*t(lsym)-t(lmin);
	tlmax = 2*t(lsym)-t(lmax);
	trmin = 2*t(rsym)-t(rmin);
	trmax = 2*t(rsym)-t(rmax);
    
	% in case symmetrized parts do not extend enough
	if tlmin(1) > t(1) | tlmax(1) > t(1)
		if lsym == indmax(1)
			lmax = fliplr(indmax(1:min(end,nbsym)));
		else
			lmin = fliplr(indmin(1:min(end,nbsym)));
		end
		if lsym == 1
			error('bug')
		end
		lsym = 1;
		tlmin = 2*t(lsym)-t(lmin);
		tlmax = 2*t(lsym)-t(lmax);
	end   
    
	if trmin(end) < t(lx) | trmax(end) < t(lx)
		if rsym == indmax(end)
			rmax = fliplr(indmax(max(end-nbsym+1,1):end));
		else
			rmin = fliplr(indmin(max(end-nbsym+1,1):end));
		end
	if rsym == lx
		error('bug')
	end
		rsym = lx;
		trmin = 2*t(rsym)-t(rmin);
		trmax = 2*t(rsym)-t(rmax);
	end 

	xlmax =x(lmax); 
	xlmin =x(lmin);
	xrmax =x(rmax); 
	xrmin =x(rmin);

	tmin = [tlmin t(indmin) trmin];
	tmax = [tlmax t(indmax) trmax];
	xmin = [xlmin x(indmin) xrmin];
	xmax = [xlmax x(indmax) xrmax];

%----------END OF boundary_conditions---------------------------------
function [instAmp,instFreq,PHI] = INST_FREQ_local(C,fs)
%%


% obtain instantaneous amplitude, frequency and phase information for IMF 'C'

ts=1/fs;

dimension=size(C);

for k=1:dimension(1)

    h=hilbert(C(k,:));

    u=real(h);
    v=imag(h);

    instAmp_temp=zeros(1,length(h));

    for n=1:length(h)   
        instAmp_temp(n)=sqrt((u(n)^2) + (v(n)^2)); 
    end
    
    instAmp(k,:)=instAmp_temp(:);
end

instAmp(1)=instAmp(2);
instAmp(end)=instAmp(end-1);




for k=1:dimension(1)
   h=hilbert(C(k,:));
   instFreq_temp=zeros(1,length(h));

   phi = unwrap(angle(h));

   PHI=(angle(h));
   for n=2:length(h)-1
       instFreq_temp(n)=(phi(n+1)-phi(n-1))/(2*ts);
   end
   
   instFreq_temp(1)=instFreq_temp(2);
   instFreq_temp(end)=instFreq_temp(end-1);
   instFreq(k,:)=instFreq_temp(:)/(2*pi);
end

instFreq=instFreq.*(double(instFreq>0));