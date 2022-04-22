function generateAllRKandLagData(data)
    [~, ~, ~, keepRsquaredValue] = validVariables();
    for j = 1:length(data.bacteria)
        allRKLagData{j} = getRKdataFromLabel(data.bacteria{j},data,keepRsquaredValue);
        allRKLagData{j}.bacterium = data.bacteria{j};
    end
    save('./data/RKandLagData.mat','allRKLagData');
end