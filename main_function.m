%     Two-step separation method to remove the topographic scatterings in dense
% array data. 
%     Matlab functions FDCTSeparate_Bayes, BalancedMultiLSAMF_new are needed to
% run this script. 
%     And the curvelet transform (Cand¨¨s et al., 2006) functions are also used, 
% which can be downloaded at http://www.curvelet.org/software.html, or you can 
% just use them presented in this package.
% 
% By Zhang et al. (2021)

clear
close all

% Two windows settings
taxis = -1:0.01:9;
Win1 = 100:600;
Win2 = 501:1001;

% parameters of the separation method
N0 = 5;
l = 70;
Norm = 2;
L = 10;
lambda1 = 0.8;
lambda2 = 1.2;
ita = 1.2;

% Filter settings
% 1(2) time-domain filtering of Real data(& Synthetic data)
% 3(4) Curvelet-domain filtering of Real data(& Synthetic data)
flagfilt = 3;   % switch
dt = 0.01;
lowf = 0.03;
highf = 4;
fs = 1/dt;
nyq = fs/2;
[B,A] = butter(4,[lowf/nyq highf/nyq]);

load('Synthetics.mat');

% If you want to plot the figure 6 in the article, you can change the ReceReal_1D to ReceReal_3D
ReceReal = ReceReal_1D(1:1001,:);
% ReceReal = ReceReal_3D(1:1001,:);
ReceSyn = ReceSyn(1:1001,:);

% weighted average
overlap = intersect(Win1,Win2);
weight1 = repmat(linspace(1,0,length(overlap))',[1 size(ReceReal,2)]);
weight2 = repmat(linspace(0,1,length(overlap))',[1 size(ReceReal,2)]);

% remove the direct P waves
Ctemp = fdct_wrapping(ReceReal,1,1);
[s1,s2] = size(Ctemp{1}{1});
Ctemp{1}{1}(1:round(200/1001*s1),:) = 0;
if flagfilt >= 3 
    Czero = fdct_wrapping(zeros(size(ReceReal)),1,1);
    nscale = length(Czero);
    for iscale = 1:nscale-1
        Czero{iscale} = Ctemp{iscale};
    end
    Ctemp = Czero;
end
tempReal = ifdct_wrapping(Ctemp,1,size(ReceReal,1),size(ReceReal,2));
Ctemp = fdct_wrapping(ReceSyn,1,1);
Ctemp{1}{1}(1:round(200/1001*s1),:) = 0;
if flagfilt == 4
    Czero = fdct_wrapping(zeros(size(ReceSyn)),1,1);
    nscale = length(Czero);
    for iscale = 1:nscale-2
        Czero{iscale} = Ctemp{iscale};
    end
    Ctemp = Czero;
end
tempSyn = ifdct_wrapping(Ctemp,1,size(ReceSyn,1),size(ReceSyn,2));

% time-domain filter
if flagfilt == 1
    for i = 1:size(tempReal,2)
        tempReal(:,i) = filtfilt(B,A,tempReal(:,i));
    end
elseif flagfilt == 2
    for i = 1:size(tempSyn,2)
        tempReal(:,i) = filtfilt(B,A,tempReal(:,i));
        tempSyn(:,i) = filtfilt(B,A,tempSyn(:,i));
    end
end

% divide into two windows
ReceReal1 = tempReal(Win1,:);
ReceReal2 = tempReal(Win2,:);
ReceSyn1 = tempSyn(Win1,:);
ReceSyn2 = tempSyn(Win2,:);
[tempd01,tempm01] = BalancedMultiLSAMF_new(ReceReal1,ReceSyn1,N0,l,Norm);
% This figure shows the whole process of the first time window
% figure;subplot(3,2,1)
% imagesc(fliplr(ReceReal1)./max(ReceReal1(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Real data')
% subplot(3,2,2)
% imagesc(fliplr(ReceSyn1)./max(ReceReal1(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('1layer Synthetic data')
% subplot(3,2,3)
% imagesc(fliplr(tempd01)./max(ReceReal1(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Filtered Result Step 1')
% subplot(3,2,4)
% imagesc(fliplr(tempm01)./max(ReceReal1(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Matched Result1 Step 1')
[d01,m01] = FDCTSeparate_Bayes(ReceReal1,tempm01,L,lambda1,lambda2,ita);
% subplot(3,2,5)
% imagesc(fliplr(d01)./max(ReceReal1(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Filtered Result Step 2')
% subplot(3,2,6)
% imagesc(fliplr(m01)./max(ReceReal1(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Matched Result1 Step 2')

[tempd02,tempm02] = BalancedMultiLSAMF_new(ReceReal2,ReceSyn2,N0,l,Norm);
% This figure shows the whole process of the second time window
% figure;subplot(3,2,1)
% imagesc(fliplr(ReceReal2)./max(ReceReal2(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Real data')
% subplot(3,2,2)
% imagesc(fliplr(ReceSyn2)./max(ReceReal2(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('1layer Synthetic data')
% subplot(3,2,3)
% imagesc(fliplr(tempd02)./max(ReceReal2(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Filtered Result Step 1')
% subplot(3,2,4)
% imagesc(fliplr(tempm02)./max(ReceReal2(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Matched Result1 Step 1')
[d02,m02] = FDCTSeparate_Bayes(ReceReal2,tempm02,L,lambda1,lambda2,ita);
% subplot(3,2,5)
% imagesc(fliplr(d02)./max(ReceReal2(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Filtered Result Step 2')
% subplot(3,2,6)
% imagesc(fliplr(m02)./max(ReceReal2(:)))
% caxis([-1 1]*5e-1)
% colormap redblue
% title('Matched Result1 Step 2')

% merge the two windows by a weighted average
d0 = zeros(size(tempReal));
m0 = zeros(size(tempSyn));
d0(1:Win1(1)-1,:) = tempReal(1:Win1(1)-1,:);
d0(Win1(1):Win2(1)-1,:) = d01(1:Win2(1)-Win1(1),:);
d0(Win1(end)+1:Win2(end),:) = d02(Win1(end)-Win2(1)+2:Win2(end)-Win2(1)+1,:);
d0(Win2(1):Win1(end),:) = d01(Win2(1)-Win1(1)+1:Win1(end)-Win1(1)+1,:) .* weight1 + d02(1:Win1(end)-Win2(1)+1,:) .* weight2;

tempd0(1:Win1(1)-1,:) = tempReal(1:Win1(1)-1,:);
tempd0(Win1(1):Win2(1)-1,:) = tempd01(1:Win2(1)-Win1(1),:);
tempd0(Win1(end)+1:Win2(end),:) = tempd02(Win1(end)-Win2(1)+2:Win2(end)-Win2(1)+1,:);
tempd0(Win2(1):Win1(end),:) = tempd01(Win2(1)-Win1(1)+1:Win1(end)-Win1(1)+1,:) .* weight1 + tempd02(1:Win1(end)-Win2(1)+1,:) .* weight2;

m0(1:Win1(1)-1,:) = tempSyn(1:Win1(1)-1,:);
m0(Win1(1):Win2(1)-1,:) = m01(1:Win2(1)-Win1(1),:);
m0(Win1(end)+1:Win2(end),:) = m02(Win1(end)-Win2(1)+2:Win2(end)-Win2(1)+1,:);
m0(Win2(1):Win1(end),:) = m01(Win2(1)-Win1(1)+1:Win1(end)-Win1(1)+1,:) .*weight1 + m02(1:Win1(end)-Win2(1)+1,:) .* weight2;

tempm0(1:Win1(1)-1,:) = tempSyn(1:Win1(1)-1,:);
tempm0(Win1(1):Win2(1)-1,:) = tempm01(1:Win2(1)-Win1(1),:);
tempm0(Win1(end)+1:Win2(end),:) = tempm02(Win1(end)-Win2(1)+2:Win2(end)-Win2(1)+1,:);
tempm0(Win2(1):Win1(end),:) = tempm01(Win2(1)-Win1(1)+1:Win1(end)-Win1(1)+1,:) .*weight1 + tempm02(1:Win1(end)-Win2(1)+1,:) .* weight2;

% Plot figure 5(6) and figure S3 in the article
traceID = 1:size(tempReal,2);
xticklabel = {};
xticklabelnull = {};
xtick = 0:10:size(tempReal,2);
for i = 1:length(xtick)
    if mod(xtick(i),50) == 0
        xticklabel{i} = num2str(xtick(i));
    else
        xticklabel{i} = '';
    end
    xticklabelnull{i} = '';
end
yticklabel = {};
ytick = -1:1:9;
for i = 1:length(ytick)
    if mod(ytick(i),2) == 0
        yticklabel{i} = num2str(ytick(i));
    else
        yticklabel{i} = '';
    end
end


subplot(4,2,1)
imagesc(traceID,taxis,fliplr(ReceReal)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('Real data')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabelnull);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
subplot(4,2,2)
imagesc(traceID,taxis,fliplr(ReceSyn)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('One-layer Synthetic data')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabelnull);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
subplot(4,2,3)
imagesc(traceID,taxis,fliplr(tempReal)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('RFs for Two-layer Model')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabelnull);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
subplot(4,2,4)
imagesc(traceID,taxis,fliplr(tempSyn)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('Predicted Scatterings')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabelnull);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
subplot(4,2,8)
imagesc(traceID,taxis,fliplr(d0)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('Filtered RFs')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabel);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
% xlabel('Seismic station');
subplot(4,2,7)
imagesc(traceID,taxis,fliplr(m0)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('Matched Scatterings')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabel);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
% xlabel('Seismic station');

subplot(4,2,6)
imagesc(traceID,taxis,fliplr(tempd0)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('Matched-filtered RFs')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabel);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
xlabel('Seismic station');

subplot(4,2,5)
imagesc(traceID,taxis,fliplr(tempm0)./max(tempReal(:)))
caxis([-1 1]*5e-1)
colormap redblue
title('Adjusted Predictions')
% set(gca,'XTick',xtick);
% set(gca,'XTicklabel',xticklabel);
ylim([-1 9])
set(gca,'YTick',ytick);
set(gca,'YTicklabel',yticklabel);
set(gca,'tickdir','out');
set(gca,'box','on');
ylabel('Time(s)');
xlabel('Seismic station');

set(gcf,'position',[0,0,700,900]);