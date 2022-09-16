function [ODdata,sugardata] = importData(fileName)

    % Import the data
    try
        [~, ~, raw] = xlsread(fileName,'Annotated');
    catch
        [~, ~, raw] = xlsread(fileName,'Sheet1');
    end
    raw = raw(:,2:49);
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells

    % Create output variable
    data = reshape([raw{:}],size(raw));
    
    sugardata = data(1,:);
    ODdata = data(2:73,:);

end