%% Local function
function climatedata = extractGdas(FilePath)
Datastore = datastore(FilePath, FileExtensions = '.txt', Type = 'tabulartext');
for f = 1 : numel(Datastore.Files)
    FilePath = Datastore.Files{f};
    % Timestamps
    DataLineStart = 7;
    DataLineEnd = 7;
    DataLineInterval = 34;
    MaxNumLines = 1e5;
    Opts = delimitedTextImportOptions(DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
        CommentStyle = ["P", ": "], ...
        ConsecutiveDelimitersRule = 'join', ...
        TrailingDelimitersRule = 'ignore');
    Cell = readcell(FilePath, Opts);
    UTC = datetime(Cell, InputFormat = 'yy MM dd HH mm');
    % Coordinate
    DataLineStart = 8;
    DataLineEnd = 8;
    DataLineInterval = 34;
    MaxNumLines = 1e5;
    Opts = delimitedTextImportOptions(DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
        CommentStyle = ["Used", "Lat: "], ...
        Delimiter = [", Lon:", " "], ...
        ConsecutiveDelimitersRule = 'join', ...
        TrailingDelimitersRule = 'ignore');
    Cell = readcell(FilePath, Opts);
    Lat = cell2mat(Cell(:, 1));
    Lon = cell2mat(Cell(:, 2));
    % 2D fields
    VariableNamesLine = 10;
    DataLineStart = 12;
    DataLineEnd = 12;
    DataLineInterval = 34;
    MaxNumLines = 1e5;
    Opts = delimitedTextImportOptions(VariableNamesLine = VariableNamesLine, ...
        DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
        Delimiter = ' ', ...
        LeadingDelimitersRule = 'ignore', ...
        ConsecutiveDelimitersRule = 'join', ...
        TrailingDelimitersRule = 'ignore');
    Field2D = readtable(FilePath, Opts);
    VariableNames = Field2D.Properties.VariableNames;
    IsExtraVariables = contains(VariableNames, 'ExtraVar');
    VariableNames = VariableNames(~IsExtraVariables);
    Field2D = Field2D(:, 2 : end);
    Field2D = varfun(@(x) str2double(x), Field2D);
    Field2D.Properties.VariableNames = VariableNames;
    Field2D = addvars(Field2D, UTC, Lat, Lon, Before = 1, NewVariableNames = {'UTC', 'LAT', 'LON'});
    % 3D fields
    VariableNamesLine = 14;
    DataLineStart = 16;
    DataLineEnd = 38;
    DataLineInterval = 34;
    MaxNumLines = 1e5;
    Opts = delimitedTextImportOptions(VariableNamesLine = VariableNamesLine, ...
        DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
        LeadingDelimitersRule = 'ignore', ...
        ConsecutiveDelimitersRule = 'join', ...
        TrailingDelimitersRule = 'ignore', ...
        Delimiter = ' ');
    Field3D = readtable(FilePath, Opts);
    VariableNames = Field3D.Properties.VariableNames;
    VariableNames = ["PRES",Field3D.Properties.VariableNames(1:end-1)];
    % IsExtraVariables = contains(VariableNames, 'ExtraVar');
    % VariableNames = VariableNames(~IsExtraVariables);
    % Field3D = Field3D(:, 2 : end);
    Field3D = varfun(@(x) str2double(x), Field3D);
    Field3D.Properties.VariableNames = VariableNames;
    % From the line below ...
    IsMissingRows = any(ismissing(Field3D), 2);
    WronglyLocatedVarNames =  {'WWND', 'RELH', 'TPOT'};
    WronglyLocatedData = Field3D{IsMissingRows, WronglyLocatedVarNames};
    MissingVarNames = {'WWND', 'RELH'};
    NonMissingVarNames = {'TPOT', 'WDIR', 'WSPD'};
    Field3D{IsMissingRows, MissingVarNames} = NaN;
    Field3D{IsMissingRows, NonMissingVarNames} = WronglyLocatedData;
    % ... until the line above is a quick & dirty workaround. This is done because WWND and RELH have missing data
    % at certain HGTS, but the readtable wrongly shifts the existing data to the left by two columns.
    % So this workaround relocates the shifted data to the correct columns. Check and validate this workaround
    % will work for all text files.
    RepeatedUTC = repelem(UTC, DataLineEnd - DataLineStart + 1);
    RepeatedLat = repelem(Lat, DataLineEnd - DataLineStart + 1);
    RepeatedLon = repelem(Lon, DataLineEnd - DataLineStart + 1);
    Field3D = addvars(Field3D, RepeatedUTC, RepeatedLat, RepeatedLon, Before = 1, NewVariableNames = {'UTC', 'LAT', 'LON'});
    % Concatenate  and sort table
    if f == 1
        Fields2D = Field2D;
        Fields3D = Field3D;
    else
        Fields2D = [Fields2D; Field2D];
        Fields3D = [Fields3D; Field3D];
    end
    if f == numel(Datastore.Files)
        Fields2D = sortrows(Fields2D, 'UTC');
        Fields3D = sortrows(Fields3D, 'UTC');
    end
end
climatedata = struct(Fields2D = Fields2D, Fields3D = Fields3D);
end