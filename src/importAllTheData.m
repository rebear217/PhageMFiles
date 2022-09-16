clear all
close all
clc

%%

load('./data/phi-secondTry-large.mat');

%%

data.igNoreSugars = [2 3];
data.sugarString = {'0','2.5','5,','10','20','40','80','160','250\mug/ml'};
data.sugarString2 = {'0\mug/ml','2.5\mug/ml','5\mug/ml','10\mug/ml','20\mug/ml','40\mug/ml','80\mug/ml','160\mug/ml','250\mug/ml'};

data.fileNames = {'2b_RYTO.xlsx','4a_RYTO.xlsx','7a_RYTO.xlsx','8a_RYTO.xlsx','9a_RYTO.xlsx',...
    '11a_RYTO.xlsx','13a_RYTO.xlsx','13b_RYTO.xlsx','14a_RYTO.xlsx',...
    '17a_RYTO.xlsx','18a_RYTO.xlsx','19a_RYTO.xlsx','21b_RYTO.xlsx',...
    '22a_RYTO.xlsx','23b_RYTO.xlsx','26a_RYTO.xlsx','27a_RYTO.xlsx','28b_RYTO.xlsx',...
    '29a_RYTO.xlsx','30a_RYTO.xlsx','41a_RYTO.xlsx','47a_RYTO.xlsx','50b_RYTO.xlsx',...
    '51a_RYTO.xlsx','51b_RYTO.xlsx',...
    '56a_RYTO.xlsx','57b_RYTO.xlsx','59a_RYTO.xlsx','60a_RYTO.xlsx',...
    '61a_RYTO.xlsx','61b_RYTO.xlsx','62b_RYTO.xlsx','66a_RYTO.xlsx',...
    '68a_RYTO.xlsx','69a_RYTO.xlsx','70a_RYTO.xlsx','70b_RYTO.xlsx',...
    '71a_RYTO.xlsx','71b_RYTO.xlsx','94a_RYTO.xlsx',...
    '95a_RYTO.xlsx','96b_RYTO.xlsx','97b_RYTO.xlsx',...
    '98a_RYTO.xlsx','99a_RYTO.xlsx','anc_RYTO.xlsx'};

data.T = 24;

data.proteinClosenessParameter = 0.005;

%%

for file = 1:length(data.fileNames)
    disp(data.fileNames{file});
    
    [data.OD{file},data.sugars{file}] = importData(['./dataRepo/rawData/OD-RYTO-runs/',data.fileNames{file}]);
    
	bacNames{file} =  data.fileNames{file}(1:3);
    
    if strcmp(bacNames{file}(3),'_')
        bacNames{file} = bacNames{file}(1:2);
	end
    if strcmp(bacNames{file}(1:2),'an')
        bacNames{file} = 'wt';
	end
end


%%

data.bacteria = bacNames;

data.infection.phi = phi;
data.infection.bacteria_names = bacteria_names;
data.infection.phage_names = phage_names;

clear file phi bacteria_names phage_names dataSet option_ug_Per_ml

%%

load('./data/allCarbonData.mat')

for k = 1:length(bacNames)
    fi = find(strcmp(bacNames{k},bacteria_names));
    data.carbonData{k} = carbonData{fi}.Cdata;
end

%%

clear bacNames bacteria_names carbonData fi k
save('./data/allData.mat');

