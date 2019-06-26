function states = segmentation_cmapss_classification(data_train,param)
%% Generate a segmentation of CMAPSS datasets into 3 states
%
% can be used for classification or clustering of time-series, as well as
% sequence prediction. The input must be in [0,1] and be globally monotonic.
%
% This file is partly related to ref [2] and to RULCLIPPER algorithm [1] available also 
% on Matlab Central. It deterministically breaks a signal into states automatically. 
% This allows to generate a "ground truth" on some signals to be used and shared in 
% classification or clustering of time-series, as well as sequence prediction
%
% Inputs: 
% data_train is a cell a time-series, each cell data_train{i} has its own length Ti
%
% param is a structure containing 
% param.part1ofsignal = 20; (default) allows to estimate the knee in the
% signal to finally get the state number 2
% param.plotfig = false; (default) % to plot or not figures
% param.thresholdResidue = 10; (default) % the threshold to decide whether
% we are in a new states
%
% Usages: 
% 1) Download RULCLIPPER algorithm (matlab central), apply the instructions in the
% readme.txt, then use:
% datasetNumber = 1; %could be 1 to 4
% [A, B, C]=health_indicators_estimation_cmapss(datasetNumber, [7 8 9 10 12 14 16 17 20 25 26], 6, true);
% This allows to get health indicators of dataset #1, see the help in 
% health_indicators_estimation_cmapss.m, B and C contains cells of data
% then 
% E=segmentation_cmapss_classification(B) or
% E=segmentation_cmapss_classification(C) will compute the segmentation into
% 3 states for health indicators estimated with the global or local method proposed in [1]
% figure,plot(E{1}) : states for instance 1
%
% 2) Second usage: You can also use it on other signals as follows:
% x=[0:0.01:5]'; E=segmentation_cmapss_classification({rand(size(x))+1.5*(1-exp(0.5*x))})
% x will be automatically normalised on [0,1] for segmentation. Not really
% validated in such manner. 
% 
% Author: E. Ramasso (emmanuel.ramasso@femto-st.fr), Dec. 2015
%
% References
%
% [1] E. Ramasso, Investigating computational geometry for failure prognostics,
% International Journal on Prognostics and Health Management, 5(5):1-18, 2014.
%
% [2] P. Cano and E. Ramasso, Ascertainment-adjusted parameter estimation approach 
% to improve robustness against misspecification of health monitoring
% methods, Mechanical Systems and Signal Processing, 2016.
%
% Version 1.1 January 9th, 7.00 am
% 

states=cell(length(data_train),1);
if nargin==1
    disp('Apply default for params')
    param.part1ofsignal = 20;
    param.plotfig = false;
    param.thresholdResidue = 10;
end
if ~isfield(param,'part1ofsignal') 
    param.part1ofsignal = 20; disp('Have change param.part1ofsignal to 20');
end    
if ~isfield( param,'plotfig') 
     param.plotfig = false; disp('Have change  param.plotfig  to false');
end    
if ~isfield(param,'thresholdResidue') 
    param.thresholdResidue  = 10; disp('Have change  param.thresholdResidue  to 10');
end    

% Go over the cell data_train
for i=1:length(data_train)
   
   x= data_train{i};
   
   % should be in [0 1] for working
   m=min(x); M=max(x); x=(x-m)./(M-m);
   
   % Take the first 20% and find a model
   r=round(param.part1ofsignal*length(x)/100);
   deb=x(1:r);
   bdeb = robustfit((1:r)'./r,deb);
   fin=x(end-r+1:end);
   bfin = robustfit((length(x)-r+1:length(x))'./length(x),fin);
   zdeb = ((1:length(x))'./length(x)) .*repmat(bdeb(2),length(x),1) + bdeb(1);
   
   if param.plotfig
       figure, plot(x)
       hold on, plot(1:length(x),zdeb,'k')
   end
   
   % apply the model for the last part, part3, this allows to find the knee
   zfin = ((1:length(x))'./length(x)) .*repmat(bfin(2),length(x),1) + bfin(1);
   
  if param.plotfig
      hold on, plot(1:length(x),zfin,'k')
      title('Beginning')
  end
  
   % residue
   r=abs(zdeb - x);
   r(1:10)=-inf;
   r = smooth(r,25);
   
   if param.plotfig
       figure,plot(r)
       title('Residue')
   end
   
   % find the time index where the residue becomes large, say around 10% of
   % the beginning
   f = find(r>=param.thresholdResidue/100);
   
   % insert states
   states{i} = 2+zeros(size(data_train{i}));
   states{i}(1:f(1))=1;
   dd = f(1);
   r=abs(zfin - x);
   r = smooth(r,25); % smooth signals in case of noise
   r(end-9:end)=+inf;
   
   if param.plotfig
       figure,plot(r)
   end
   
   f = find(r<=param.thresholdResidue/100);
   f(find(f<=dd+param.thresholdResidue/100*length(r)))=[];
   cont=isempty(f);
   states{i}(f(1):end)=3;

   if param.plotfig
       figure,plot(states{i})
       pause
       close all
   end
   
   assert(length(unique(states{i}))==3)
end

figure,  
for i=1:length(data_train), 
    r=rand(1,3); 
    subplot(211), hold on, plot(data_train{i},'color',r)
    subplot(212), hold on, plot(states{i},'color',r), 
end
subplot(211), xlabel('Time unit'), ylabel('Health indicators'), title('Data'), axis tight
subplot(212), ylabel('States'), xlabel('Time unit'), axis tight

%figure_pdf_cropped(gcf,'s1','png')

