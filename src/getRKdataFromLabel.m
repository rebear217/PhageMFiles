function RKdata = getRKdataFromLabel(label,data,keepRsquaredValue,parallelFlag)

    if nargin < 4
        parallelFlag = 1;
    end
    if nargin < 3 || isempty(keepRsquaredValue)
        keepRsquaredValue = 0.7;
    end
    
    warning('off')
    
    if parallelFlag
        try
            nc = feature('numcores');
            if matlabpool('size') == 0
                matlabpool('open',min(12,2*nc));
            end
        catch
            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(poolobj)
                poolsize = 0;
                disp('creating parallel pool');
                parpool;
            else
                poolsize = poolobj.NumWorkers;
            end
        end
    end
    
    sugars = unique(data.sugars{1});
    sugars = sugars(~isnan(sugars));
    
    Rsquareddata = [];
    Rdata = [];
    growthRateData = [];
    Lagdata = [];
    Lagdata2 = [];
    Lagdata3 = [];
    Kdata = [];
    Rsedata = [];
    Ksedata = [];
    growthRateSEData = [];
    allSugars = [];

    parfor j = 1:length(sugars)
        FitData{j} = logisticFitToODfromLabel(sugars(j),data,label,0);
    end
    
	for j = 1:length(sugars)
        Rsquareddata = [Rsquareddata FitData{j}.Rsquared];
        Rdata = [Rdata FitData{j}.R];
        growthRateData = [growthRateData FitData{j}.growthRate];
        Kdata = [Kdata FitData{j}.K];
        Rsedata = [Rsedata FitData{j}.Rse];
        Ksedata = [Ksedata FitData{j}.Kse];
        growthRateSEData = [growthRateSEData FitData{j}.growthRateSE];
        
        % FitData{j}.lagTimeMeasure1 is not being used
        
        %model-derived lag:
        Lagdata = [Lagdata FitData{j}.lagTimeMeasure2];

        %this is the F-test-derived lag time from seeking non-constant regression to OD:
        Lagdata2 = [Lagdata2 FitData{j}.lagTimeMeasure3];
       
        %this is NAN (a place-holder for more lag stats):
        Lagdata3 = [Lagdata3 FitData{j}.lagTimeMeasure4];
        allSugars = [allSugars sugars(j)*ones(size(FitData{j}.R))];
    end

    badFits = length(find(Rsquareddata < keepRsquaredValue));
    Fits = length(Rsquareddata);
    disp([label,' has ',num2str(badFits),' / ',num2str(Fits),' rejected logistic OD datafits']);
    
    keep = (Rsquareddata > keepRsquaredValue) .* ~isnan(Rsedata) .* ~isnan(Ksedata) .* ~isnan(Rdata) .* ~isnan(Kdata) .* ~isnan(growthRateData);
    
    Rsquareddata = Rsquareddata(keep>0);
    Rdata = Rdata(keep>0);
    Kdata = Kdata(keep>0);
    Rsedata = Rsedata(keep>0);
    Ksedata = Ksedata(keep>0);
    allSugars = allSugars(keep>0);
    Lagdata = Lagdata(keep>0);
    Lagdata2 = Lagdata2(keep>0);
    Lagdata3 = Lagdata3(keep>0);
    growthRateData = growthRateData(keep>0);
    growthRateSEData = growthRateSEData(keep>0);
    
    Yielddata = Kdata./allSugars;
    Yieldsedata = Ksedata./allSugars;
    Yielddata = Yielddata(~isnan(Yielddata));
    Yieldsedata = Yieldsedata(~isnan(Yielddata));
   
    RKdata.sugars = sugars;
    RKdata.Rsquareddata = Rsquareddata;
    RKdata.Rdata = Rdata;
    RKdata.growthRateData = growthRateData;
    RKdata.Kdata = Kdata;
    RKdata.Rsedata = Rsedata;
    RKdata.Ksedata = Ksedata;
    RKdata.growthRateSEData = growthRateSEData;
    RKdata.allSugars = allSugars;
    RKdata.Yielddata = Yielddata;
    RKdata.Yieldsedata = Yieldsedata;
    RKdata.Lagdata = Lagdata;
    RKdata.Lagdata2 = Lagdata2;
    RKdata.Lagdata3 = Lagdata3;
end