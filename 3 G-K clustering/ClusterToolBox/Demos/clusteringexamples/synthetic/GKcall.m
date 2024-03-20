clear all
close all
path(path,'..\..\..\FUZZCLUST')
%data set
data.X = xlsread('C:\Users\WangLu\Desktop\Code\3 G-K clustering\para_example.xlsx');

%parameters
param.c=4;
param.m=2;
param.e=1e-6;
param.ro=ones(1,param.c);
param.val=1;
%normalization
data=clust_normalize(data,'range');

result = GKclust(data,param);
result = validity(result,data,param);

plot(data.X(:,1),data.X(:,2),'b.',result.cluster.v(:,1),result.cluster.v(:,2),'ro');
hold on
%draw contour-map
new.X=data.X;
eval=clusteval(new,result,param);
result.validity