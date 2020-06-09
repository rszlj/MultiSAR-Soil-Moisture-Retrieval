function [mv] = soilTimeSeries(meanmv,stdmv)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
len=length(meanmv);
mv=meanmv+randn(1,len).*stdmv;
mv(mv>0.44)=0.44;
mv(mv<0.04)=0.04;
[mv]=sort(mv,'descend');
end

