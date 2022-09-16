function I = getLocationLabels(name,nameSet)
    I = find(cellfun(@(s)strcmpi(s,name),nameSet));
end