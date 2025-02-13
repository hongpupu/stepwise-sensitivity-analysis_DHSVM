clear all
close all
path(path,'..\..\..\FUZZCLUST')

%the data
data.X = xlsread('I:\��һ��\01 ����\Clustering\test.xlsx');

%parameters
param.c=3;
param.m=2;
param.e=1e-6;
param.ro=ones(1,param.c);
%normalization
data=clust_normalize(data,'range');

result = FCMclust(data,param);
param.val=1;
result = validity(result,data,param);
plot(data.X(:,1),data.X(:,2),'b.',result.cluster.v(:,1),result.cluster.v(:,2),'ro');
hold on
%draw contour-map
new.X=data.X;
eval=clusteval(new,result,param);
result.validity