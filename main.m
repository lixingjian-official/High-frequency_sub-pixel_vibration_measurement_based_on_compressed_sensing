clear
%% Add MatlabPyrTools and Natsortfiles
addpath(fullfile('matlabPyrTools'));addpath('natsortfiles\');

%% Read Image Files
currentDirectory = pwd;
dataDir = uigetdir();
pics = dir(strcat(dataDir,'\','*.tif'));
[~,ind] = natsort({pics.name});
pics = pics(ind);
n = length(pics)-1;
S_orgin = getmotionsignal(dataDir,1:n);

%% Images Sampled
Compression_Rate = 0.3;
m = floor(Compression_Rate*n);
index = randperm(n,m);
Psi = inv(fft(eye(n,n)));     
A = Psi(index,:); 

%% Recover Vibration Signals from Sampled Images
S = getmotionsignal(dataDir,index);
f2 = S.aligned_horizontal; % Calculate the Horizontal Vibration Signal
f2 = S.aligned_vertical; % Calculate the Vertical Vibration Signal

%% Recover Vibration Signals by Compression Sensing
cvx_begin quiet % Convex Optimization Toolbox,CVX
    variable x(n) complex;
    minimize( norm(x,1) );
    subject to
      A*x == f2 ;
      norm( x, Inf ) <= 1e3
cvx_end
['CVX finished !']

x3=ADMM_L1_reconstruct(1e3*A,1e3*f2,100,10,100); % Alternating Direction Method of Multipliers,ADMM
['ADMM finished !']

S_orgin.alignedout = fft(S_orgin.aligned_horizontal); % Ground Truth of Horizontal Vibration Signal
S_orgin.alignedout = fft(S_orgin.aligned_vertical); % Ground Truth of Vertical Vibration Signal
x4 = S_orgin.alignedout;

%% Plot Signals in Frequency Dimension
end_num = n/2-1;num = 0:end_num;
figure('Units','centimeter','Position',[10 14 18 6]);
plot(num(2:end),abs(x(2:end_num+1))/max(abs(x(2:end_num+1))),'color',[0 0.4470 0.7410]);hold on 
plot(num(2:end),abs(x3(2:end_num+1))/max(abs(x3(2:end_num+1))),'color',[0.8500 0.3250 0.0980]);hold on 
plot(num(2:end),abs(S_orgin.alignedout(2:end_num+1))/max(abs(S_orgin.alignedout(2:end_num+1))), ...
    'color',[0.9290 0.6940 0.1250]);hold on 
plot([1 499], [0.1 0.1],'--k');
yticks([0 0.1 0.5 1]);
xlabel('Frequency(Hz)');ylabel('Amplitude (Normalized)')
set(gca,'FontName','Times New Roman');
legend('CVX','ADMM','GT');
box off     
ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  
set(ax1,'XTick', [],'YTick', []);   
hold off

%% Plot Signals in Time Dimension
lengthout = 100;num = 1:lengthout;
CVXout = real(ifft(x));ADMMout = real(ifft(x3));GTout = real(ifft(S_orgin.alignedout));

figure('Units','centimeter','Position',[10 5 18 6]);
plot(num,mapminmax(CVXout(1:lengthout)',0,1), 'color',[0 0.4470 0.7410]);hold on 
plot(num,mapminmax(ADMMout(1:lengthout)',0,1), 'color',[0.8500 0.3250 0.0980]);hold on 
plot(num,mapminmax(GTout(1:lengthout)',0,1),'color',[0.9290 0.6940 0.1250])
yticks([0 0.5 1]); 
xlabel('Time(ms)');ylabel('Amplitude (Normalized)')
set(gca,'FontName','Times New Roman');
legend('CVX','ADMM','GT');
box off     
ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', []);   
hold off

%% Calculate MSE
CVXout = mapminmax(abs(ifft(x))',0,1);
ADMMout = mapminmax(abs(ifft(x3))',0,1);
GTout = mapminmax(abs(ifft(x4))',0,1);
outmse1 = mse(GTout-CVXout)
outmse3 = mse(GTout-ADMMout)