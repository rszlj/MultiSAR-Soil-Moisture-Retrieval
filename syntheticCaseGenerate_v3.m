%Synthetic data sets (full simulation in frequency and incidence angle)
clear;
load("C:\Users\zljnj\Dropbox\Matlab_Codes\EnsembleLearning_soil moisture\Cubes4BareLC\forwardCube_IIEM");
load("C:\Users\zljnj\Dropbox\Matlab_Codes\EnsembleLearning_soil moisture\SyntheticInput\OzNet_sm.mat"); % the soil moiustre of OzNet
caseNum=200; % a total of 200 cases were simulated
inputCase=cell(caseNum,1);
mvGT=cell(caseNum,1);
%Roughness
sig=rand(200,1)*2.5+0.5; %RMS height in [0.5, 3] cm
cl=rand(200,1)*30+5; %correlation length [5,35] cm
%Radar Configuration
fre=[1.26,5.4]; % frequency in GHz
incid=[23,35]; % incidence angle in deg

% generate the synthetic dataset of case A
for kkk=1:caseNum
rng(kkk,'twister'); % set random number generator
%soil moisture
meanmv=nanmean(OzNet_sm)/100;
stdmv=nanstd(OzNet_sm)/100;
[mv] = soilTimeSeries(meanmv,stdmv); % generate the simulated time series of soil moisture according to the average and standard deviation of OzNet observations
%Forward sigmaNull
dateNum=length(meanmv); % the number of dates, or the number of time instances
tID4L=1:dateNum;
tbs=nan(dateNum,3);
%L band@23 deg
for i=1:dateNum
[ tbs(i,1),tbs(i,2),tbs(i,3) ] =bsFromCubes( sig(kkk),cl(kkk),mv(i),incid(1),forwardCube_IIEM.cubeL ); % interpolate the sigmanull from LUTs
end
lia=ones(dateNum*3,1)*incid(1); % incidence angle
sideID=ones(dateNum*3,1)*1; % radar look direction 1 left, 2 right
fr=ones(dateNum*3,1)*1; % freqnecy id, 1: L-band, 2:C-band, 3: X-band
bsID=repmat([1,2,3],dateNum,1);
bsID=bsID(:); % polarizaiton id, 1: HH, 2:HV, 3:VV
tID=repmat(tID4L',3,1);
tID=tID(:);% date sequence
temp1=[tbs(:),lia,bsID,tID,sideID,fr]; % 96 (3 pol * 8 dates * 2 frequency * 2 incidence angle)*6(sigmanull, incidence angle, polarization, sequence, look direction, frequency)


%L band@35 deg
for i=1:dateNum
[ tbs(i,1),tbs(i,2),tbs(i,3) ] =bsFromCubes( sig(kkk), cl(kkk), mv(i),incid(2),forwardCube_IIEM.cubeL );
end
lia=ones(dateNum*3,1)*incid(2);
sideID=ones(dateNum*3,1)*2;
fr=ones(dateNum*3,1)*1;
bsID=repmat([1,2,3],dateNum,1);
bsID=bsID(:);
tID=repmat(tID4L',3,1);
tID=tID(:);
temp2=[tbs(:),lia,bsID,tID,sideID,fr];

%C band@23 deg
for i=1:dateNum
[ tbs(i,1),tbs(i,2),tbs(i,3) ] =bsFromCubes(sig(kkk),cl(kkk),mv(i),incid(1),forwardCube_IIEM.cubeC );
end
lia=ones(dateNum*3,1)*incid(1);
sideID=ones(dateNum*3,1)*1;
fr=ones(dateNum*3,1)*2;
bsID=repmat([1,2,3],dateNum,1);
bsID=bsID(:);
tID=repmat(tID4L',3,1);
tID=tID(:);
temp3=[tbs(:),lia,bsID,tID,sideID,fr];

%C band@35 deg
for i=1:dateNum
% [ tbs(i,1),tbs(i,2),tbs(i,3) ] =bsFromCubesBatch_Bag( sig,mv(i),incid(2),5.41 );
[ tbs(i,1),tbs(i,2),tbs(i,3) ] =bsFromCubes( sig(kkk),cl(kkk),mv(i),incid(2),forwardCube_IIEM.cubeC );

end
lia=ones(dateNum*3,1)*incid(2);
sideID=ones(dateNum*3,1)*2;
fr=ones(dateNum*3,1)*2;
bsID=repmat([1,2,3],dateNum,1);
bsID=bsID(:);
tID=repmat(tID4L',3,1);
tID=tID(:);
temp4=[tbs(:),lia,bsID,tID,sideID,fr];

bs=[temp1;temp2;temp3;temp4]; %each simulation have four parts L-band at 23 and 35 deg and C-band ata 23 and 35 deg
% case A: only contain speckle noise
Nlevel=0.7;
[m,n]=size(bs(:,1));
bs(:,1)=bs(:,1)+normrnd(0,0.7,m,n);
inputCase{kkk}=bs;
mvGT{kkk}=mv;
end
save('CaseA','inputCase','sig','cl', 'mvGT')% save the sigmanull, rms height, correlation length, soil moisture

% case B: adding 0.5 to HH and -0.5 to HV
clear
load('CaseA.mat')
for i=1:200
    temp=inputCase{i};
    temp(temp(:,3)==1,1)=temp(temp(:,3)==1,1)+0.5;
    temp(temp(:,3)==2,1)=temp(temp(:,3)==2,1)-0.5;
    inputCase{i}=temp;
end
save('CaseB','inputCase','sig','cl', 'mvGT')% save the sigmanull, rms height, correlation length, soil moisture

% Case C: adding 5, -1, -1.5 and -2 dB to C-band 23°, C-band 35°, L-band 23° and L-band 35°, respectively
clear
load('CaseA.mat')
for i=1:200
    temp=inputCase{i};
    index=temp(:,2)*10+temp(:,6);
    temp(index==231,1)=temp(index==231,1)-2;
    temp(index==232,1)=temp(index==232,1)-1.5;
    temp(index==351,1)=temp(index==351,1)-1;
    temp(index==352,1)=temp(index==352,1)-5;
    inputCase{i}=temp;
end
save('CaseC','inputCase','sig','cl', 'mvGT')% save the sigmanull, rms height, correlation length, soil moisture
