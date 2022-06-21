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