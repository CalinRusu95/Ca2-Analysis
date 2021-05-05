%Clear console, figures and workspace
clear all;
close all;

%Open and prepare results file
res=fopen('C:\Users\calin\Google Drive\Proiecte\Calciu\CaImag\CaPeak\results.txt', 'w+');
specs = fopen('C:\Users\calin\Google Drive\Proiecte\Calciu\CaImag\CaPeak\spectre.txt', 'w+');
%fprintf(res,'%9s %10s %5s %15s %15s %13s %16s %19s %15s %11s %11s %13s \n','#File', 'Duration', 'Type', 'Asymmetry', 'IncreaseSpeed1', 'IncreaseSpeed2', 'DecreaseSpeed1', 'DecreaseSpeed2', 'Area', 'Amplitude1', 'Amplitude2', 'PeaksDistance');

%Import input files and input parameters
InputPath = 'C:\Users\calin\Google Drive\Proiecte\Calciu\CaImag\CaPeak\Input\2021.04.23\14.04.2021-1, 6 Gy, CF, FF, puls atp\48 h\6 Gy FF\';

cd(InputPath);
files = dir('**');
files(1:2) = [];
FilesNum = numel(files);
cplot = ceil(sqrt(FilesNum));
InSigs = [];
corrRes = [];
%opts = setvaropts(opts,'HexType','uint64');

for i =1:FilesNum %For all input files
    
    %Get dataset size and parameters
    FileAdrs=strcat(InputPath,files(i).name);
    name = files(i).name;
    opts = detectImportOptions(FileAdrs);
    %opts = setvartype(opts,{'Var1'},'int64');
    %opts = setvartype(opts,{'Var2'},'double');
    %opts.HexType = {'uint64'};
    A = readtable(FileAdrs,opts);
    [M,N] = size(A);
    X=A(:,1); %Time
    spectra=A(:,2:N); %Intensities
        
    %Set spectra limits and starting time
    fracini=0.13;
    %fracini=0.01;
    
    fracfin=0.1;
    %fracfin=0.01;
    
    pkfrac=0.05;
    %pkfrac=0.05;
    
    promfrac=0.1;
    %promfrac=0.1;
    
    MPD=5000;
    %MPD=3000;
    
    tstart=120000; %first transient
    %tstart=905000; %3 min 30' second transient
    %tstart=965000; %4 min 30' second transient
    %tstart=1302000; %10 min second transient
    
    %Compute spectra features
    for j=1:N-1
        
        %Smooth the spectra and get gradient and maximum value
        spectra = table2array(spectra);
        Yout = smooth(spectra,11,'sgolay');
        spectra(:,j)=Yout;
        intd=gradient(Yout);
        [c,l]=max(Yout);
        
        %Fit data with a Gauss distribution
        [m,s] = normfit(Yout);
        gaussian = normpdf(Yout,m,s);
        
        %Get maximum time, value and increase speed
        X = table2array(X);
        tmax(j)=X(l);
        peakmax(j)=c;
        vmax(j)=max(intd);
        vmin(j)=min(intd);
        
        first=1; %Check if first
        last=3; %Check if last
        
        %Determine initial and final points (time, value) of the peak
        for k=1:M-1
            if Yout(k+1)>fracini*c && first==1    
                tini(j)=(X(k)+X(k+1))/2;
                kini=k;
                first=2;
            end
            
            if X(k)>tmax(j)
                if Yout(k+1)<fracfin*c &&  last==3
                    tfin(j)=(X(k)+X(k+1))/2;
                    kfin=k;
                    last=4;
                end;
            end 
        end;
        
        xpeak=X(kini:kfin); %Peak X axis
        peak=Yout(kini:kfin); %Peak Y axis
        xpeak2 = typecast(xpeak,'double');
        aria(j)=trapz(xpeak2, peak);
    end;
    
    [pks,locs,widths,proms] = findpeaks(peak,xpeak,'MinPeakHeight',pkfrac*c,'MinPeakDistance',MPD,'MinPeakProminence',promfrac*c);
  
    %Plot all input signals in separate files
    figure
	plot(X, spectra)
    title(name)
    
    %Check the spectra and peak for each input file
    figure(i)
    subplot(1,2,1)
    plot(X, spectra)
    hold on
    plot(xpeak,peak)
    title(strcat(name,' - signal and peak'))
    subplot(1,2,2)
    plot(xpeak, peak, 'r')
    findpeaks(peak,xpeak,'MinPeakHeight',pkfrac*c,'MinPeakDistance',MPD,'MinPeakProminence',promfrac*c)
    title(strcat(name,' - peak'))
    
    %Plot all input signals
    figure(FilesNum+2)
    subplot(cplot,cplot,i)
    plot(X, spectra)
    title(name)
    
    %Plot all peaks
    figure(FilesNum+3)
    subplot(cplot,cplot,i)
    plot(X, spectra)
    hold on
    plot(xpeak,peak)
    title(name)
    hold off
    
    %Compute increase and decay time
    Tincrease=tmax-tini;
    Tdekay=tfin-tmax;
    
    %Compute results
    Ttotal=tfin-tini;
    Latency=tini-tstart;
    Asymmetry=Tdekay./Tincrease;
    Increase1Speed=vmax;
    Decrease1Speed=abs(vmin);
    Increase2Speed=0;
    Decrease2Speed=0;
    PeakArea=aria;
    Peak1Amplitude=peakmax;
    Peak2Amplitude=0;
    PeakNumber = length(pks);
    Type = 1;
    PeaksDistance = 0;
    
    if PeakNumber==2 && pks(1)>pks(2)
        Type = 2;
        PeaksDistance = abs(locs(1)-locs(2));
        Peak1Amplitude = pks(1);
        Peak2Amplitude = pks(2);
        gradis=sort(intd);
        Increase2Speed=gradis(end-1);
        Decrease2Speed=abs(intd(2));
    else
        if PeakNumber==2 && pks(1)<pks(2)
            Type = 3;
            PeaksDistance = abs(locs(1)-locs(2));
            Peak1Amplitude = pks(1);
            Peak2Amplitude = pks(2); 
            gradis=sort(intd);
            Increase2Speed=gradis(end-1);
            Decrease2Speed=abs(intd(2));
        end
    end
    
    param = [Type, Ttotal, Latency, Asymmetry, Increase1Speed, Increase2Speed, Decrease1Speed, Decrease2Speed, PeakArea, Peak1Amplitude Peak2Amplitude PeaksDistance];
    corrRes(i,1:7)= [Ttotal, Latency, Asymmetry, Increase1Speed, Decrease1Speed, PeakArea, Peak1Amplitude];
    
    %Print results to file
    fprintf(res, '%9s %3.0f %10.1f %10.1f %15.6f %15.6f %15.6f %15.6f %15.6f %12.4f %11.4f %11.4f %13.0f \n', name, param);
    
    figure(FilesNum+4)
    [H,AX,BigAx,P,PAx] = plotmatrix(corrRes);
    set(gca,'xlim',[0 5])
    ylabel(AX(1,1),'Ttotal','FontSize',10,'Rotation',90)
    ylabel(AX(2,1),'Latency','FontSize',10,'Rotation',90)
    ylabel(AX(3,1),'Asymmetry','FontSize',10,'Rotation',90)
    ylabel(AX(4,1),'Inc. Speed','FontSize',10,'Rotation',90)
    ylabel(AX(5,1),'Dec. Speed','FontSize',10,'Rotation',90)
    ylabel(AX(6,1),'Area','FontSize',10,'Rotation',90)
    ylabel(AX(7,1),'Amplitude','FontSize',10,'Rotation',90)
    
    xlabel(AX(7,1),'Ttotal','FontSize',10,'Rotation',0)
    xlabel(AX(7,2),'Latency','FontSize',10,'Rotation',0)
    xlabel(AX(7,3),'Asymmetry','FontSize',10,'Rotation',0)
    xlabel(AX(7,4),'Inc. Speed','FontSize',10,'Rotation',0)
    xlabel(AX(7,5),'Dec. Speed','FontSize',10,'Rotation',0)
    xlabel(AX(7,6),'Area','FontSize',10,'Rotation',0)
    xlabel(AX(7,7),'Amplitude','FontSize',10,'Rotation',0)
    %[R,p]=corr(paramSub1, 'type', 'Pearson')
    %[R,p]=corr(param, 'type', 'Spearman');
    %[R,p]=corr(param, 'type', 'Kendall');
    
end

%Close results file
fclose(res);
fclose(specs);