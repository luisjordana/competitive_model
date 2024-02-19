
function [imx, hilph]=halfcyclehilbert(x,pk,nmlz)

% Calculates imaginary part of hilbert transform of signal x for each half
% cycle (imx), and according Hilbert phase (hilph). This procedure takes into account that 
% a signal does not always nicely oscillate around zero for every half cycle and that 
% amplitudes may differ widely over half cycles. It is particularly fit for higly variable
% signals.
% inputs: x= oscillating signal; pk = peak indices determined with e.g. findpeaks.m ;
% nmlz = optional normalization (default == 0; advise: just keep it zero)
% 



if nargin==2
    nmlz=0;
end

if size(x,2)>1
    x=x';  % x should be column array
end

% create arrays with NaNs:
imx=ones(length(x),1)*NaN;
xtemp=imx;
n=10;% n repetitions of reconstructed halfcycle signal

for i=1:length(pk)-1;
    %create half-cycle signal
    tempsig=x(pk(i):pk(i+1));
    % signal centre to zero:
    tempsig=tempsig-mean(tempsig);
    % optional: normalize tempsig to amplitude = 1   (haalt niet heel erg veel anders uit)
    if nmlz==1
        tempsig=tempsig-mean([tempsig(1) tempsig(end)]);
        tempsig=tempsig/range(tempsig);
    end
    xtemp(pk(i):pk(i+1))=tempsig;
    
    % create longer signal for proper hilbert transform:
    temp=[tempsig; flipud(tempsig(2:end-1))];
    temp= repmat(temp',1,n);
    % do hilbert transform and take imaginary part:
    imtemp=imag(hilbert(temp));
    % cut centre piece out:
    imx(pk(i):pk(i+1)-1)=imtemp((length(tempsig)-1)*n+1:(length(tempsig)-1)*(n+1));
    
    clear temp tempsig imtemp
end
hilph=cart2pol(xtemp,imx);