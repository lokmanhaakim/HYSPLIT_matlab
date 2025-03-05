function Data = HYSPLITRun(Location,Lt,Ln,Time,Period,Map,currentfile,savedfiles,Met_type)

arguments
    Location (1,:) string
    Lt (1,:) double
    Ln (1,:) double
    Time (1,:) datetime
    Period (1,:) double
    Map table
    currentfile (1,1) string
    savedfiles string = ""
    Met_type string = "gdas1"
end



Location = arrayfun(@(x) erase(x," "),Location);

% Set all directory C:\Users\Admin\Desktop\nuclear-wind
workdir = fullfile("C:","HYSPLIT","working");
execdir = fullfile("C:","HYSPLIT","exec");
metdir = fullfile("C:","Users","Admin","Desktop","MetData");%C:\Users\Admin\Desktop\MetData
CompletedSimulation = 0;

if ~strcmp(savedfiles,"")

    savedir = fullfile("C:","Users","Admin","Desktop","nuclear-wind",savedfiles);

    if ~exist(savedir, 'dir')
        mkdir(savedfiles);
    end

end

Outertic = tic;
for l = 1:numel(Location)

    Location(l)
    Data.(Location(l)) = [];

    for st =1:numel(Time)
        Location (l)
        Time(st)
        Innertic = tic;

        tic
        try
            web = sprintf("https://api.open-meteo.com/v1/elevation?latitude=%.2f&longitude=%.2f",Lt(l),Ln(l));
            Altitude = webread(web);
            Altitude = Altitude.elevation;

        catch
            Altitude = 10;
        end

        if Altitude > 100
            Altitude = 50;
        end

        %[8,27,50,44.1,10,10,11]; 11.6739° N, 108.                 %["NortheastMonsoon" "SouthwestMonsoon"]; %NE = 22 dec 2023 , SW = 19 May 2023 , Inter-monsoon = 1 Apr 2023

        end_time = Time(st)+hours(Period);
        diff_week = week(end_time)-week(Time(st));

        if Met_type == "NCEP"
            Period_interval= Time(st):calmonths(1):end_time+1;

        elseif Met_type == "gfs025"
            Period_interval= Time(st):days(1):end_time+1;

        else
            Met_type = "gdas1";
            Period_interval= Time(st):days(7):end_time+1;
        end

        emission = 3.030e16;
        emssion_rate=1.0;
        hrs = "06";
        min = "00";

        % Climate data
        cd(workdir);
        climatedata.Fields3D = [];
        climatedata.Fields2D = [];

        metfileName = "";

        for s = 1:numel(Period_interval)
            met_dir = fullfile(metdir,Met_type,"/");

            if Met_type == "NCEP"
                Mnth = sprintf("%02d",month(Time(st)));
                metfileName(s,:) = "RP"+string(year(Time(st))+Mnth+".gbl");

            elseif Met_type == "gfs025"
                metfileName(s,:) = sprintf("%d%02d%02d_gfs0p25",year(Period_interval(s)),month(Period_interval(s)),day(Period_interval(s)));

            else
                metfileName(s,:) = "gdas1."+lower(string(month(Period_interval(s),"shortname")))+string(mod(year(Period_interval(s)),100))+".w"+string(ceil(day(Period_interval(s))./7));

            end


            % Define the command with the current coordinates
            %C:\\Users\\USER\\OneDrive\\Desktop\\Master\\nuclear-wind\\met\\NARR
            command = sprintf('"%s\\profile.exe" -d%s\\ -f\\%s  -y%.2f -x%.2f -t1 -w1',execdir,met_dir,metfileName(s),Lt(l), Ln(l));
            % Run the command
            system(command);
            if ans == 900 % 900 means it has error
                disp(metfileName(s)+" not exist")
                disp("The file is not exist, please wait for " + metfileName(s) +" downloading process...");
                ftpobj = ftp("ftp://arlftp.arlhq.noaa.gov");

                if Met_type == "NCEP"
                    cd(ftpobj,"archives/reanalysis");
                elseif Met_type == "gfs025"
                    cd(ftpobj,"archives/gfs0p25");
                else
                    cd(ftpobj,"archives/gdas1");
                end

                mget(ftpobj,metfileName(s),met_dir);
                command = sprintf('"C:\\HYSPLIT\\exec\\profile.exe" -dC:\\Users\\USER\\OneDrive\\Desktop\\Master\\nuclear-wind\\met\\NARR\\%s\\ -f\\%s  -y%.2f -x%.2f -t1 -w1',metfileName,Lt(l), Ln(l));
                % Run the command
                system(command);
            end

            if Met_type == "NCEP"
                climatedataT = readNCEP('profile.txt');
            elseif Met_type == "gfs025"
                climatedataT = readGFS('profile.txt');
            else
                climatedataT = readGDAS('profile.txt');
            end


            climatedata.Fields3D = sortrows([climatedata.Fields3D;climatedataT.Fields3D],"UTC");
            climatedata.Fields2D = sortrows([climatedata.Fields2D;climatedataT.Fields2D],"UTC");

        end

        % HYSPLIT - concentration run (still not solved)

        % hysplitPath = 'C:\HYSPLIT\';
        % system("C:\HYSPLIT\exec\hycs_std.exe")
        % hycs_std
        %"C:\HYSPLIT\working\CONTROL" -- > untuk parameter setup
        ctFile = './CONTROL';
        output_name = "tdump";
        output_dir = "./";

        % met_dir = fullfile(metdir,string(year(Period_interval(s))))+"\";
        fid = fopen(ctFile, 'w');

        Yr = mod(year(Time(st)),100); %last 2 digit
        Yr = string(Yr);
        Mnth = sprintf("%02d",month(Time(st)));
        Dy=sprintf("%02d",day(Time(st)));
        Hr = sprintf("%02d",hour(Time(st)));
        num_loc = numel(Lt(l));
        mix_layer = "1000.0";
        top = "10000.0";

        period_accident = Period;
        vertical_layer = 0;

        num_polutant = 1;
        nuclide = "C137";
        emission = string(emission);
        emssion_rate= string(emssion_rate);
        start_hourAccident = "17";

        fprintf(fid, '%s %s %s %s\n',Yr,Mnth,Dy,Hr);
        fprintf(fid, '%d\n', num_loc);
        fprintf(fid, '%0.1f %0.1f %0.1f\n',Lt(l),Ln(l),Altitude);
        fprintf(fid,'%d\n',period_accident);
        fprintf(fid,'%d\n',vertical_layer);
        fprintf(fid,'%s\n',top);
        fprintf(fid,'%d\n',numel(metfileName));

        for m = 1:numel(metfileName)
            fprintf(fid,'%s\n',met_dir);
            fprintf(fid,'%s\n',metfileName(m));
        end

        fprintf(fid,'%s\n',output_dir);
        fprintf(fid,'%s\n',output_name);
        fclose(fid);

        command = fullfile(execdir,'hyts_std.exe');                       %"C:\HYSPLIT\exec\hyts_std.exe";
        system(command)%
        command = "C:\HYSPLIT\exec\trajplot.exe -a1 -iC:\HYSPLIT\working\tdump";
        system(command)
        Traj = readtable(fullfile(workdir,"GIS_traj_ps_01.txt"));
        Traj.Properties.VariableNames = ["Pressure","Longitude","Latitude","Height"];
        Traj(end,:) = [];
        Traj.Distance = deg2km(distance('gc',Lt(l),Ln(l),Traj.Latitude,Traj.Longitude));
        Traj.Bearing = azimuth(Lt(l),Ln(l),Traj.Latitude,Traj.Longitude);

        ctFile = './CONTROL';
        % met_dir = sprintf("C:\\Users\\USER\\OneDrive\\Desktop\\Master\\nuclear-wind\\met\\NARR\\%s",string(year(Period_interval(s))));
        fid = fopen(ctFile, 'w');

        Yr = mod(year(Time(st)),100); %last 2 digit
        Yr = string(Yr);
        Mnth = sprintf("%02d",month(Time(st)));
        Dy=sprintf("%02d",day(Time(st)));
        Hr = sprintf("%02d",hour(Time(st)));
        num_loc = numel(Lt(l));
        mix_layer = "1000.0";
        top = "10000.0";

        period_accident = Period;
        vertical_layer = 0;

        num_polutant = 1;
        nuclide = "C137";
        emission = string(emission);
        emssion_rate= string(emssion_rate);
        start_hourAccident = "17";

        num_grids = "1";
        centre_lat = "0.0";
        centre_lon = "0.0";
        spacing_lat = "0.05";
        spacing_lon = "0.05";
        span_lat = "50.0";
        span_lon = "50.0";
        output_dir = "./";
        output_name = "cdump";
        vert_levels = "4";
        height_down = "0";
        height_up = ["10" "20" "30" "40"];
        type_calc = "00";
        hrs = "06";
        min = "00";

        deposition ="1";
        particlediameter = "5.0 6.0 1.0";
        constant1 = "0.006 0.0 0.0 0.0 0.0";
        constant2 = "0.0 8.0E-05 8.0E-05";
        decay = "10960.0";
        re_factor = "0.0";

        stime = datetime(year(Time(st)),month(Time(st)),day(Time(st)),"Format","uu MM dd");
        stimeaccident = datetime(year(Time(st)),month(Time(st)),day(Time(st)),"Format","uuuu MM dd");
        samplestart = datetime(year(Time(st)),month(Time(st)),day(Time(st)),"Format","uuuu MM dd");

        fprintf(fid, '%s %s %s %s\n',Yr,Mnth,Dy,Hr);
        fprintf(fid, '%d\n', num_loc);
        fprintf(fid, '%0.1f %0.1f %0.1f\n',Lt(l),Ln(l),Altitude);
        fprintf(fid,'%d\n',period_accident);
        fprintf(fid,'%d\n',vertical_layer);
        fprintf(fid,'%s\n',top);
        fprintf(fid,'%d\n',numel(metfileName));

        for m = 1:numel(metfileName)
            fprintf(fid,'%s\n',met_dir);
            fprintf(fid,'%s\n',metfileName(m));
        end

        fprintf(fid,'%d\n',num_polutant);
        fprintf(fid,'%s\n',nuclide);
        fprintf(fid,'%s\n',emission);
        fprintf(fid,'%s\n',emssion_rate);
        fprintf(fid, '%s %s %s %s\n',Yr,Mnth,Dy,Hr);
        fprintf(fid,'%s\n',num_grids);
        fprintf(fid,'%s ',centre_lat);
        fprintf(fid,'%s\n',centre_lon);
        fprintf(fid,'%s ',spacing_lat);
        fprintf(fid,'%s\n',spacing_lon);
        fprintf(fid,'%s ',span_lat);
        fprintf(fid,'%s\n',span_lon);
        fprintf(fid,'%s\n',output_dir);
        fprintf(fid,'%s\n',output_name);
        fprintf(fid,'%s\n',vert_levels);
        fprintf(fid,'%s ',height_down);
        fprintf(fid,'%s %s %s %s\n',height_up(1),height_up(2),height_up(3),height_up(4));
        fprintf(fid,"00 00 00 00 00\n");
        fprintf(fid,"00 00 00 00 00\n");
        fprintf(fid,'%s %s %s \n',type_calc,hrs,min);
        fprintf(fid,'%s\n',deposition);
        fprintf(fid,'%s\n',particlediameter);
        fprintf(fid,'%s\n',constant1);
        fprintf(fid,'%s\n',constant2);
        fprintf(fid,'%s\n',decay);
        fprintf(fid,'%s\n',re_factor);
        fclose(fid);

        command = "C:\HYSPLIT\exec\hycs_std.exe";
        system(command)

        command = "C:\HYSPLIT\exec\con2asc.exe -d, -iC:\HYSPLIT\working\cdump -ocdump -s1";
        system(command)
        %
        Alt = Altitude(1);
        opts = delimitedTextImportOptions(Delimiter=",",LeadingDelimitersRule="ignore",VariableNamesLine=1);

        Conc = readtable(fullfile(workdir,"cdump.txt"),opts); %insert text concentration
        Conc(1,:) = [];
        Conc = varfun(@(x) str2double(x),Conc);
        Conc.UTC = datetime(Conc.Fun_YEAR,Conc.Fun_MO,Conc.Fun_DA,Conc.Fun_HR,0,0);
        Conc(:,["Fun_DA" "Fun_HR" "Fun_MO" "Fun_YEAR"]) = [];
        Csidx = find(contains(Conc.Properties.VariableNames,"Fun_C13700"));
        strCs = string(Conc.Properties.VariableNames(Csidx));

        for strCs = strCs
            Conc.(strCs)(Conc.(strCs)==0,:) = NaN;
        end

        isEmpty = isnan(Conc.Fun_C13700000) & isnan(Conc.Fun_C13700010) & isnan(Conc.Fun_C13700020) & isnan(Conc.Fun_C13700030);
        Conc(isEmpty,:) = [];

        Conc.Distance = deg2km(distance('gc',Lt(l),Ln(l),Conc.Fun_LAT,Conc.Fun_LON)); %l s
        Conc.Bearing = azimuth(Lt(l),Ln(l),Conc.Fun_LAT,Conc.Fun_LON);
        Conc = movevars(Conc,"UTC","Before","Fun_LAT");

        for sh = 1:height(Map)
            % Extract the polygon vertices for the current region
            polygonLon = Map.Lon{sh};
            polygonLat = Map.Lat{sh};

            % Check if points are inside the polygon for all Conc points
            Region = inpolygon(Conc.Fun_LON, Conc.Fun_LAT, polygonLon, polygonLat);

            % Update Conc.Region based on the current region (sh)
            Conc.Region(Region) = sh;
        end
        %
        % Data.(Location(l)).StartDate(st) = Time(st);
        % Data.(Location(l)).ClimateData{st} = climatedata;
        % Data.(Location(l)).Concentration{st} = Conc;
        % Data.(Location(l)).Trajectory{st} = Traj;

        Variable = Conc.Properties.VariableNames(4:end-3);
        ldata = groupsummary(Conc,"UTC",["max","mean","std"],Variable);
        u = unique(Conc.UTC);
        [ConcCell,RegionCell,TrajCell,CentroidCell] = deal(cell(numel(u), 1));
        for d = 1:numel(u)
            ConcCell{d,:} = Conc(Conc.UTC==u(d),:);
            TrajCell{d,:} = Traj;
            RegionCell{d,:} = unique(ConcCell{d}.Region);
            RegionCell{d,:} = nonzeros(unique(cat(2,ConcCell{d}.Region)));
            [c1,c2] = centroid(polyshape(ConcCell{d}.Fun_LAT,ConcCell{d}.Fun_LON)); %,Trajectory = Traj
            CentroidCell{d,:} = [c1,c2];
        end
        ldata.Concentration = ConcCell;
        ldata.RegionCell = RegionCell;
        ldata.Trajectory = TrajCell;
        ldata.Centroid = CentroidCell;

        % Data.(Location(l)).ldata{st} = ldata;

        % Store data
        Data.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(st)),month(Time(st)),day(Time(st)),hour(Time(st)))).climatedata = climatedata;
        Data.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(st)),month(Time(st)),day(Time(st)),hour(Time(st)))).Concentration = Conc;
        Data.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(st)),month(Time(st)),day(Time(st)),hour(Time(st)))).Trajectory = Traj;
        Data.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(st)),month(Time(st)),day(Time(st)),hour(Time(st)))).ldata = ldata;

        if ~strcmp(savedfiles,"")

            fileName = sprintf("%s.mat",Location(l));
            SavefilePath = fullfile(savedir,fileName);

            save(SavefilePath,"-struct","Data",Location(l),"-v7.3")

        end

        CompletedSimulation = CompletedSimulation+1;
        disp("Percentage completed :" + (CompletedSimulation*100)./(numel(Time)*(numel(Location)))+"%")
        Innertoc = toc(Innertic)

    end


end

Outertoc = toc(Outertic)

cd(currentfile)

    function climatedata = readGDAS(FilePath)
        % Timestamps
        DataLineStart = 6;
        DataLineEnd = 6;
        DataLineInterval = 34;
        MaxNumLines = 1e5;
        Opts = delimitedTextImportOptions(DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
            CommentStyle = ["P", ": "], ...
            ConsecutiveDelimitersRule = 'join', ...
            TrailingDelimitersRule = 'ignore');
        Cell = readcell(FilePath, Opts);
        UTC = datetime(Cell, InputFormat = 'yy MM dd HH mm');
        % Coordinate
        DataLineStart = 7;
        DataLineEnd = 7;
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
        VariableNamesLine = 9;
        DataLineStart = 11;
        DataLineEnd = 11;
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
        VariableNamesLine = 13;
        DataLineStart = 15;
        DataLineEnd = 37;
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

        climatedata = struct(Fields2D = Field2D, Fields3D = Field3D);

    end

    function climatedata = readGFS(FilePath)

        % Timestamps
        DataLineStart = 6;
        DataLineEnd = 6;
        DataLineInterval =66;
        MaxNumLines = 1e5;
        Opts = delimitedTextImportOptions(DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
            CommentStyle = ["P", ": "], ...
            ConsecutiveDelimitersRule = 'join', ...
            TrailingDelimitersRule = 'ignore');
        Cell = readcell(FilePath, Opts);
        UTC = datetime(Cell, InputFormat = 'yy MM dd HH mm');
        % % Coordinate
        DataLineStart = 7;
        DataLineEnd = 7;
        DataLineInterval = 66;
        MaxNumLines = 1e5;
        Opts = delimitedTextImportOptions(DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
            CommentStyle = ["Used", "Lat: "], ...
            Delimiter = [", Lon:", " "], ...
            ConsecutiveDelimitersRule = 'join', ...
            TrailingDelimitersRule = 'ignore');
        Cell = readcell(FilePath, Opts);
        Lat = cell2mat(Cell(:, 1));
        Lon = cell2mat(Cell(:, 2));
        % % % 2D fields
        VariableNamesLine = 9;
        DataLineStart = 11;
        DataLineEnd = 11;
        DataLineInterval = 66;
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
        % % % 3D fields
        VariableNamesLine = 13;
        DataLineStart = 15;
        DataLineEnd = 28;
        DataLineInterval = 66;
        MaxNumLines = 1e5;
        Opts = delimitedTextImportOptions(VariableNamesLine = VariableNamesLine, ...
            DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
            LeadingDelimitersRule = 'ignore', ...
            ConsecutiveDelimitersRule = 'join', ...
            TrailingDelimitersRule = 'ignore', ...
            Delimiter = ' ');
        Field3D = readtable(FilePath, Opts);
        VariableNames = Field3D.Properties.VariableNames;
        Field3D(:,7) = [];
        VariableNames(strcmp(VariableNames,'PRES')) = [];
        VariableNames = ["PRES",VariableNames(1:end-1)];
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

        climatedata = struct(Fields2D = Field2D, Fields3D = Field3D);

    end

    function climatedata = readNCEP(FilePath)

        % Timestamps
        DataLineStart = 6;
        DataLineEnd = 6;
        DataLineInterval =28;
        MaxNumLines = 1e5;
        Opts = delimitedTextImportOptions(DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
            CommentStyle = ["P", ": "], ...
            ConsecutiveDelimitersRule = 'join', ...
            TrailingDelimitersRule = 'ignore');
        Cell = readcell(FilePath, Opts);
        UTC = datetime(Cell, InputFormat = 'yy MM dd HH mm');
        % Coordinate
        DataLineStart = 7;
        DataLineEnd = 7;
        DataLineInterval = 28;
        MaxNumLines = 1e5;
        Opts = delimitedTextImportOptions(DataLines = [(DataLineStart : DataLineInterval : MaxNumLines)', (DataLineEnd : DataLineInterval : MaxNumLines)'], ...
            CommentStyle = ["Used", "Lat: "], ...
            Delimiter = [", Lon:", " "], ...
            ConsecutiveDelimitersRule = 'join', ...
            TrailingDelimitersRule = 'ignore');
        Cell = readcell(FilePath, Opts);
        Lat = cell2mat(Cell(:, 1));
        Lon = cell2mat(Cell(:, 2));
        % % 2D fields
        VariableNamesLine = 9;
        DataLineStart = 11;
        DataLineEnd = 11;
        DataLineInterval = 28;
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
        % % 3D fields
        VariableNamesLine = 13;
        DataLineStart = 15;
        DataLineEnd = 28;
        DataLineInterval = 28;
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

        climatedata = struct(Fields2D = Field2D, Fields3D = Field3D);

    end

end