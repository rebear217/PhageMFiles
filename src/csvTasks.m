function csvTasks()

    % run from root folder:
    % first find the list of files to process to CSV
    
    d = dir('./dataRepo/**');
    dN = length(d);
    fileList = {};
    count = 1;
    for j = 1:dN
        thisFile = d(j).name;
        if (length(thisFile) > 4)
            ending = thisFile(end-3:end);
            if strcmp(ending,'.xls') || strcmp(ending,'xlsx')
                f = [d(j).folder,'/',thisFile];
                fileList{count} = f;
                count = count + 1;
            end
        end
    end
    
    % Now process them:
    
    fN = length(fileList);
    for j = 1:fN
        disp(['Completed ',num2str(j),' of ',num2str(fN),' CSV files.'])
        fName = fileList{j};
        sheets = sheetnames(fName);
        lS = length(sheets);
        for k = 1:lS
            sheet = sheets{k};
            Table = readtable(fName,'Sheet',sheet);

            sheetName = removePunc(sheet);
            header = split(fName,'.');
            header = header{1};

            if lS > 1
                sheetFolder = [header,'-csvSheets/'];
                if ~exist(sheetFolder,'dir')
                    mkdir(sheetFolder)
                end
                fileName = [sheetFolder,sheetName,'.csv'];
            else
                fileName = [header,'.csv'];
            end
            writetable(Table,fileName)
        end
    end

end

function result = removePunc(str)
    keeperIndexes = str == ' ' | (str>='0' & str<='9') | (str>='a' & str<='z') | (str>='A' & str<='Z');
    result = str(keeperIndexes);
end
