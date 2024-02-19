function [realpeaks,amplitudes,minpeaks,maxpeaks]=peakfind(x,minampl,pct)

% Finds real peaks in oscillating signal x
% depending on amplitude
%
% Minampl is the minimal amplitude of  
% (positive) peak(i) and (negative) peak(i+1)
% Is the amplitude below this value, than the cycle
% is not considered a real movement cycle
%
% if 'pct' = 1, 'minampl' will be considered as a percentage
%     of the mean amplitude
%
% Choose value of minampl yourself by looking at 
% data: it also depends on the movement frequency (not
% included in this script)
%
% Harjo de Poel, 8-7-2004
% - adaptation 'first ptemp value', Harjo de Poel, 10-5-2013

if nargin==2
    pct=0;
end

dx=diff(x);
nd=length(dx);
countermin=0;
countermax=0;


% for i=2:nd-1
%    if ((dx(i)>0 & dx(i-1)<0)|(dx(i)==0 & (dx(i-1)<0 & dx(i+1)>0)) | ...
%            (dx(i)<0 & dx(i-1)>0)|(dx(i)==0 & (dx(i-1)>0 & dx(i+1)<0))) & ...
%            abs(x(i)-x(i-1))>minampl
%        counter=counter+1;
%        peaks(counter,1)=i; %= t of peaks
%    end
% end



for i=2:nd-1
   if (dx(i)>0 & dx(i-1)<0)|(dx(i)==0 & (dx(i-1)<0 & dx(i+1)>0))
       countermin=countermin+1;
       minima(countermin,1)=i; %=tmin
   end
   if (dx(i)<0 & dx(i-1)>0)|(dx(i)==0 & (dx(i-1)>0 & dx(i+1)<0))
       countermax=countermax+1;
       maxima(countermax,1)=i; %=tmax
   end
end

%minima(diff(minima)==0)=[];
%maxima(diff(maxima)==0)=[];

ptemp=sort([maxima(:,1); minima(:,1)]);
ampstemp=abs(diff(x(ptemp)));

% peaks=peakstemp;
%meanamp=mean(abs(x(ptemp(1:end-1))-x(ptemp(2:end))));
meanamp=mean(ampstemp);
if pct==1
    minampl=minampl*meanamp;
end

% first make sure the first peak is one that precedes an oscillation of size > minapl
fa=find(ampstemp>minampl);
ptemp=ptemp(fa(1):end);
peaks(1)=ptemp(1);

for j=2:length(ptemp)
    amp=abs(x(ptemp(j-1))-x(ptemp(j:end)));
    ok=find(amp>minampl);
    if ~isempty(ok)
        peaks(j)=ptemp(j-1+ok(1));
        if ~isempty (find(peaks(1:j-1)==peaks(j)))
            peaks(j)=NaN;
        end
    else
        peaks(j)=NaN;
    end
end
f=find(isnan(peaks));
peaks(f)=[];
peaks=sort(peaks);
amps=x(peaks);
peaks1=peaks;
for k=1:length(peaks)-1
    if ~isempty (find(maxima==peaks(k)))
        if ~isempty (find(maxima==peaks(k+1)))
            [ampmax,maxpeaki]=max([amps(k) amps(k+1)]);
            peaks(k-maxpeaki+2)=NaN;            
        end
    elseif ~isempty (find(minima==peaks(k)))
        if ~isempty (find(minima==peaks(k+1)))
            [ampmin,minpeaki]=min([amps(k) amps(k+1)]);
            peaks(k-minpeaki+2)=NaN;
        end
    end
end
f2=find(isnan(peaks));
peaks(f2)=[];

realpeaks=peaks;
amplitudes=x(peaks);
maxpeaks=realpeaks(ismember(realpeaks,maxima));
minpeaks=realpeaks(ismember(realpeaks,minima));
