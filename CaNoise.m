%Clear console, figures and workspace
clear all;
close all;

%Open and prepare results file
res=fopen('\path_to\results.txt', 'w+');
specs = fopen('\path_to\spectre.txt', 'w+');
%fprintf(res,'%9s %10s %5s %15s %15s %13s %16s %19s %15s %11s %11s %13s \n','#File', 'Duration', 'Type', 'Asymmetry', 'IncreaseSpeed1', 'IncreaseSpeed2', 'DecreaseSpeed1', 'DecreaseSpeed2', 'Area', 'Amplitude1', 'Amplitude2', 'PeaksDistance');

%Import input files and input parameters
InputPath = '\path_to\input_data_folder\';
cd(InputPath);
files = dir('**');
files(1:2) = [];
FilesNum = numel(files);
cplot = ceil(sqrt(FilesNum));
InSigs = [];
corrRes = [];

for i =1:FilesNum %For all input files
%for i =1:1 %For 2 input files
    
    %Get dataset size and parameters
    FileAdrs=strcat(InputPath,files(i).name);
    name = files(i).name;
    A = importdata(FileAdrs);
    [N,M] = size(A);
    t=A(:,1); %Time
    dt = mean(diff(t)); 
    fs = 1 / dt; % Hz
    t = t/fs;    
    x=A(:,2:M); %Intensities
        
    %% detrend the signal
    % organize question dialog menu about the detrending
%     quest = 'Do you want to detrend the signal?';
%     dlgtitle = 'Detrending';
%     btn1 = 'Yes, detrend the signal';
%     btn2 = 'No, do not detrend the signal';
%     defbtn = btn1;
%     answer = questdlg(quest, dlgtitle, btn1, btn2, defbtn);
%     % normalize the signal
% 
%     switch answer
%         case btn1
%         % detrend the signal    
%         x = detrend(x);                             
%         case btn2
%         % do not detrend the signal
%     end
    
    %normalize the signal
%     %organize question dialog menu about the normalization
%     quest = 'What type of normalization do you want?';
%     dlgtitle = 'Normalization';
%     btn1 = 'Normalize the signal to unity peak';
%     btn2 = 'Normalize the signal to unity RMS-value';
%     btn3 = 'Do not normalize the signal';
%     defbtn = btn1;
%     answer = questdlg(quest, dlgtitle, btn1, btn2, btn3, defbtn);
%     %normalize the signal
%     
%     switch answer
%         case btn1
%         %normalize to unity peak
%         x = x/max(abs(x));
%         case btn2
%         %normalize to unity RMS-value
%         x = x/std(x); 
%         case btn3
%         %do not normalize the signal
%     end
    
    %% Plot signal oscilogram
    figure(1)
    subplot(3,2,1)
    plot(t, x, 'r')
    xlim([0 max(t)])
    ylim([-1.1*max(abs(x)) 1.1*max(abs(x))])
    grid minor
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8)
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    title('Oscillogram of the noise signal') 

    %% Plot signal periodogram
    % calculate noise PSD using Welch's periodogram
    winlen = round(N/100);
    win = blackman(winlen, 'periodic');
    hop = round(winlen/4);
    nfft = round(2*winlen);
    [PSD, f] = pwelch(x, win, winlen-hop, nfft, fs, 'onesided', 'psd');
    PSD = 10*log10(PSD);
    %figure(2)
    subplot(3,2,2)
    semilogx(f, PSD, 'r', 'LineWidth', 1.5)
    grid minor
    %xlim([100 10000])
    %ylim([-120 0])
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8)
    xlabel('Frequency (Hz)')
    ylabel({'Magnitude (I^{2}/Hz)'})
    title('Power Spectral Density of the noise signal')
    
    %% Plot signal spectrogram using time-frequency analysis
    [~, F, T, STPSD] = spectrogram(x, win, winlen-hop, nfft, fs, 'psd');
    STPSD = 10*log10(STPSD);
    %figure(3)
    subplot(3,2,3)
    surf(T, F, STPSD)
    shading interp
    axis tight
    box on
    view(0, 90)
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Spectrogram of the signal')
    [~, cmax] = caxis;
    caxis([max(-120, cmax-90), cmax])
    hClbr = colorbar;
    set(hClbr, 'FontName', 'Times New Roman', 'FontSize', 8)
    ylabel(hClbr, 'Magnitude (I^{2}/Hz)')

    %% Plot signal histogram
    %figure(4)
    subplot(3,2,4)
    hHist = histogram(x, round(sqrt(N/10)), 'FaceColor', 'r');
    xlim([-1.1*max(abs(x)) 1.1*max(abs(x))])
    ylim([0 1.1*max(get(hHist, 'Values'))])
    grid minor
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8)
    xlabel('Amplitude (V)')
    ylabel('Number of samples')
    title('Histogram of the noise signal')

    %% Autocorrelation function estimation
    [Rx, lags] = xcorr(x, 'coeff');
    tau = lags/fs;
    %figure(5)
    subplot(3,2,[5,6])
    plot(tau, Rx, 'r')
    grid minor
    xlim([-max(tau) max(tau)])
    ylim([1.1*min(Rx), 1.1])
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8)
    xlabel('Delay (s)')
    ylabel('Autocorrelation coefficient')
    title('Correlogram of the noise signal')
    line([-max(abs(tau)) max(abs(tau))], [0.05 0.05], 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--')
    legend('noise correlogram', '5 % level')

%     %% signal stationarity estimation test
%     if all(isstationary(x))
%         disp('According the stationarity estimation test the time series is stationary.')
%     else
%         disp('According the stationarity estimation test the time series is non-stationary.')
%     end
% 
    %% Statistics
    
    % minimum and maximum values
    maxval = max(x);
    minval = min(x);
%     disp(['Max value = ' num2str(maxval)])
%     disp(['Min value = ' num2str(minval)])
 
    % DC and RMS values
    u = mean(x);
    s = std(x);
    skew = skewness(x);
    kurt = skewness(x);
    disp(['Mean value = ' num2str(u)])
    disp(['RMS value = ' num2str(s)])
    
    % Skewness and kurtosis
    disp(['Skewness value = ' num2str(skewness(x))])
    disp(['Kurtosis value = ' num2str(kurtosis(x))])
    
    % Dynamic range
    DR = 20*log10(max(abs(x))/min(abs(nonzeros(x))));
    disp(['Dynamic range DR = ' num2str(DR) ' a.u.'])
   
    % Crest factor
    CF = 20*log10(max(abs(x))/s);
    disp(['Crest factor CF = ' num2str(CF) ' a.u.'])
    
    % Autocorrelation time
    ind = find(Rx>0.05, 1, 'last');
    RT = (ind-N)/fs;
    %disp(['Signal duration = ' num2str(max(t)) ' s'])
    disp(['Autocorrelation time = ' num2str(RT) ' ms'])
    
    disp('---------------------------------------------')
    
    param = [maxval, minval, u, s, skew, kurt, DR, CF, RT];
    % Maxim  Minim  Medie  SD  Skewness  Kurtosis  Interval dinamic  Factor de creasta  Timp de autocorelare
    format = '%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n';
    fprintf(res, format, param);
    
    figname = erase(files(i).name, '.csv');
    savefolder = strcat('\path_to\plot_output_folder\', figname,'.png');
    saveas(figure(1),[savefolder]);
   
    commandwindow
    
end

%Close results file
fclose(res);
fclose(specs);
