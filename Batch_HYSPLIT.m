%% Setup locations
clc;clearvars;
NPPStation = readtable("C:\Users\Admin\Desktop\nuclear-wind\NPPStationPRIS.xlsx");
NPPStation = NPPStation(~isnan(NPPStation.Longitude),:);
[~,idxa] = unique(NPPStation(:,1:2));
NPPStation = NPPStation(idxa,:);
Location = NPPStation.Location;
Location = arrayfun(@(x) regexprep(x, '[^a-zA-Z]', ''),Location);
Location = string(Location);
Lat = NPPStation.Latitude;
Lon = NPPStation.Longitude;
Time = [datetime(2021,3,14),datetime(2022,3,14)];
Altitude = zeros(length(Lon),1);        

for l = 1:numel(Lat)

    try
        web = sprintf("https://api.open-meteo.com/v1/elevation?latitude=%.2f&longitude=%.2f",Lat(l),Lon(l));
        AltitudeWeb = webread(web);
        Altitude(l) = AltitudeWeb.elevation;

    catch
        Altitude(l) = 10;
    end

end


% Period = 24*7;
NewData = [];
%% Setup emission rate (Terada et al.,2020)
EmissionRate = readtable("EmissionRateTime.xlsx","Sheet","Table015 (Page 6)","VariableNamingRule","preserve");
EmissionRate.Properties.VariableNames(1) = "UTC";
StartDate = Time(1);
Duration = EmissionRate.("Duration(h)"); 
Emission1 =  EmissionRate.("137Cs");
Emission2 = EmissionRate.("131I");

UTCExplosion = StartDate + hours(cumsum(Duration));
UTCExplosion = [StartDate;UTCExplosion(1:end-1)];

EmissionTable = table(UTCExplosion,Duration,Emission1,Emission2,'VariableNames',["UTC","Duration","137Cs","131I"]);
EmissionTable.YYYY = year(EmissionTable.UTC);
EmissionTable.MM = month(EmissionTable.UTC);
EmissionTable.DD = day(EmissionTable.UTC);
EmissionTable.HH = hour(EmissionTable.UTC);
RunTimes = hours(EmissionTable.UTC(end)-EmissionTable.UTC(1:end-1));
RunTimes = [RunTimes;0];

EndDate = EmissionTable.UTC(end);

cd("C:\HYSPLIT\working");
EmissionTable = EmissionTable(:,["YYYY" "MM" "DD" "HH" "137Cs" "131I"]);
Emission1 = arrayfun(@(x) sprintf("%.1E",x),Emission1);
Emission2 = arrayfun(@(x) sprintf("%.1E",x),Emission2);
EmissionTable.Properties.VariableNames(end-1:end) = ["Cs-137","I-131p"];
EmissionTable(:,end-1) = table(Emission1);
EmissionTable(:,end) = table(Emission2);

EmissionTable=EmissionTable(:,["YYYY" "MM" "DD" "HH"]);
EmissionTable = arrayfun(@(x) sprintf("%02d",x),EmissionTable{:,:});
EmissionTable = [EmissionTable,RunTimes,Emission1,Emission2];
writematrix(EmissionTable,"Allday.txt","Delimiter"," ");
EmissionTable(2:end+1,:) = EmissionTable;
EmissionTable(:,5)=[];
EmissionTable(1,:) = ["YYYY" "MM" "DD" "HH" "Cs-137" "I-131p"];
writematrix(EmissionTable,"TeradaEmission.txt","Delimiter"," ");
%% Met directory 

Period = RunTimes(1);
metdir = fullfile("C:","Users","Admin","Desktop","MetData");

%% Batch file 
%  metfileName(s,:) = "gdas1."+lower(string(month(Period_interval(s),"shortname")))+string(mod(year(Period_interval(s)),100))+".w"+string(ceil(day(Period_interval(s))./7));
cd("C:\HYSPLIT\working")
ctFile = './Helloworld.bat';
CompletedSimulation = 0;

for l = 1:numel(Location)
    for t = 1:numel(Time)
        Innertic = tic;
        delete('SG*');
        delete('CG*');
        delete('DG*');
        delete('TG*');

        temp_met = [];

        StartDate = Time(t);
        Duration = EmissionRate.("Duration(h)");
        Emission1 =  EmissionRate.("137Cs");
        Emission2 = EmissionRate.("131I");

        UTCExplosion = StartDate + hours(cumsum(Duration));
        UTCExplosion = [StartDate;UTCExplosion(1:end-1)];

        EmissionTable = table(UTCExplosion,Duration,Emission1,Emission2,'VariableNames',["UTC","Duration","137Cs","131I"]);
        EmissionTable.YYYY = year(EmissionTable.UTC);
        EmissionTable.MM = month(EmissionTable.UTC);
        EmissionTable.DD = day(EmissionTable.UTC);
        EmissionTable.HH = hour(EmissionTable.UTC);
        RunTimes = hours(EmissionTable.UTC(end)-EmissionTable.UTC(1:end-1));
        RunTimes = [RunTimes;0];

        EndDate = EmissionTable.UTC(end);

        cd("C:\HYSPLIT\working");
        EmissionTable = EmissionTable(:,["YYYY" "MM" "DD" "HH" "137Cs" "131I"]);
        Emission1 = arrayfun(@(x) sprintf("%.1E",x),Emission1);
        Emission2 = arrayfun(@(x) sprintf("%.1E",x),Emission2);
        EmissionTable.Properties.VariableNames(end-1:end) = ["Cs-137","I-131p"];
        EmissionTable(:,end-1) = table(Emission1);
        EmissionTable(:,end) = table(Emission2);

        EmissionTable=EmissionTable(:,["YYYY" "MM" "DD" "HH"]);
        EmissionTable = arrayfun(@(x) sprintf("%02d",x),EmissionTable{:,:});
        EmissionTable = [EmissionTable,RunTimes,Emission1,Emission2];
        writematrix(EmissionTable,"Allday.txt","Delimiter"," ");
        EmissionTable(2:end+1,:) = EmissionTable;
        EmissionTable(:,5)=[];
        EmissionTable(1,:) = ["YYYY" "MM" "DD" "HH" "Cs-137" "I-131p"];
        writematrix(EmissionTable,"TeradaEmission.txt","Delimiter"," ");

        str_Start = extractAfter(sprintf("%d%02d%02d%02d",year(StartDate),month(StartDate),day(StartDate),hour(StartDate)),2);
        str_End = extractAfter(sprintf("%d%02d%02d%02d",year(EndDate),month(EndDate),day(EndDate),hour(EndDate)),2);

        end_time = Time(t)+hours(Period);
        diff_week = week(end_time)-week(Time(t));
        Period_interval= Time(t):days(7):end_time+1;
        for s = 1:numel(Period_interval)
            met_dir = fullfile(metdir,"gdas1","/");
            metfileName(s,:) = "gdas1."+lower(string(month(Period_interval(s),"shortname")))+string(mod(year(Period_interval(s)),100))+".w"+string(ceil(day(Period_interval(s))./7));
            temp_met = [temp_met;met_dir;metfileName(s,:)];
        end

        for tM = 1:numel(temp_met)
            str_temp_met(tM,:) = sprintf("echo %s                      >>CONTROL",temp_met(tM));
        end

        fid = fopen(ctFile, 'w');

        if fid == -1
            error('Failed to open batch file.');
        end

        batchCommands = ["@echo off";...
            "setLocal EnableDelayedExpansion";...
            "set DIR=c:";...
            "set PGM=%DIR%\hysplit";...
            "   ";...
            "echo -90.0   -180.0  lat/lon of lower left corner   >ASCDATA.CFG";...
            "echo 1.0     1.0     lat/lon spacing in degrees    >>ASCDATA.CFG";...
            "echo 180     360     lat/lon number of data points >>ASCDATA.CFG";...
            "echo 2               default land use category     >>ASCDATA.CFG";...
            "echo 0.2             default roughness length (m)  >>ASCDATA.CFG";...
            "echo '%PGM%\bdyfiles\'  directory of files         >>ASCDATA.CFG"
            " ";...
            "echo ^&SETUP               >SETUP.CFG";...
            "echo maxdim = 2,          >>SETUP.CFG";...
            "echo delt = 5.0,          >>SETUP.CFG";...
            "echo khmax = 24,          >>SETUP.CFG";...
            "echo numpar = 1500,      >>SETUP.CFG";...
            "echo maxpar = 10000000,      >>SETUP.CFG";...
            "echo /                    >>SETUP.CFG";...
            "  ";...
            "REM ------------------------------------------------";...
            "   ";...
            "for /F ""tokens=1,2,3,4,5,6,7"" %%A in (Allday.txt) do (";...
            "    set YEAR=%%A";...
            "    set MONTH=%%B";...
            "    set DAY=%%C"
            "    set HOUR=%%D";...
            "    set RUN=%%E";...
            "    set SA=%%F";...
            "    set SB=%%G";...
            "echo !DAY!";...
            "echo !HOUR!";...
            "echo !SA!";...
            "echo !SB!";...
            sprintf("    echo %s %02d !DAY! !HOUR!            >CONTROL",extractBefore(str_Start,3),month(StartDate));...
            "    echo 2                       >>CONTROL";...
            sprintf("    echo %.2f %.2f %.2f    >>CONTROL",Lat(l),Lon(l),Altitude(l));...
            sprintf("    echo %.2f %.2f %.2f    >>CONTROL",Lat(l),Lon(l),Altitude(l));...
            "echo !RUN!"
            "    echo !RUN!                   >>CONTROL";...
            "    echo 0                       >>CONTROL";...
            "    echo 10000.0                 >>CONTROL";...
            sprintf("    echo %d                       >>CONTROL",numel(metfileName));...
            sprintf("    %s\n",str_temp_met);
            "    echo 2                       >>CONTROL";...
            "    echo RNUC                    >>CONTROL";...
            "    echo 1.0                     >>CONTROL";...
            "    echo !RUN!                     >>CONTROL";...
            "    echo 00 00 00 00 00          >>CONTROL";...
            "    echo NGAS                    >>CONTROL";...
            "    echo 1.0                     >>CONTROL";...
            "    echo !RUN!                    >>CONTROL";...
            "    echo 00 00 00 00 00          >>CONTROL";...
            "    echo 1                       >>CONTROL";...
            sprintf(("    echo %.2f %.2f              >>CONTROL"),Lat(l),Lon(l));...
            "    echo 0.05 0.05               >>CONTROL";...
            "    echo 20.0 30.0               >>CONTROL";...
            "    echo ./                      >>CONTROL";...
            "    echo TG_03!DAY!!HOUR!             >>CONTROL";...
            "    echo 2                       >>CONTROL";...
            "    echo 0 100                   >>CONTROL";...
            "    echo 00 00 00 00 00          >>CONTROL";...
            "    echo 00 00 00 00 00          >>CONTROL";...
            "    echo 00 03 00                >>CONTROL";...
            "    echo 2                       >>CONTROL";...
            "    echo 1.0 1.0 1.0             >>CONTROL";...
            "    echo 0.001 0.0 0.0 0.0 0.0   >>CONTROL";...
            "    echo 0.0 8.0E-05 8.0E-05     >>CONTROL";...
            "    echo 0.0                     >>CONTROL";...
            "    echo 0.0                     >>CONTROL";...
            "    echo 0.0 0.0 0.0             >>CONTROL";...
            "    echo 0.0 0.0 0.0 0.0 0.0     >>CONTROL";...
            "    echo 0.0 0.0 0.0             >>CONTROL";...
            "    echo 0.0                     >>CONTROL";...
            "    echo 0.0                     >>CONTROL";...
            "  ";...
            "    %PGM%\exec\hycs_std";...
            "   ";...
            "    REM ---- terminate day hour loop";...
            ")";...
            ")";...
            " ";...
            " ";...
            "rem ---------------------------------------------";...
            "  ";...
            "%PGM%\exec\condecay -1:1:11025.8:C137 -2:1:8.02330:I131 +eTeradaEmission +oCG_ +t031205";...
            " ";...
            "for /F ""tokens=1,2,3,4,5,6,7"" %%A in (Allday.txt) do (";...
            "set DAY=%%C";...
            "set HOUR=%%D";...
            "  %PGM%\exec\con2asc -d, -iCG_03!DAY!!HOUR! -oCG_03!DAY!!HOUR! -s1"
            "  %PGM%\exec\con2rem -iCG_03!DAY!!HOUR! -oDG_03!DAY!!HOUR! -s1 -t0 -d1 -aactivity_unit.txt";...
            "  %PGM%\exec\concsum -iDG_03!DAY!!HOUR! -oSG_03!DAY!!HOUR! -l -pDOSE";...
            "  echo Finished: !DAY! !HOUR!";...
            ")";...
            ")";...
            "  ";...
            "dir /b SG_??????.bin >merg_list.txt";...
            "%PGM%\exec\conmerge -imerg_list.txt -ofdnpp_total.bin";...
            "    ";...
            "%PGM%\exec\conavgpd -ifdnpp_total.bin -ofdnpp_sums.bin "+sprintf("-a%s -b%s -r1",str_Start,str_End);... %adjust str_StartDate,str_EndDate
            "    ";...
            "%PGM%\exec\con2asc -d, -ifdnpp_total.bin -ofdnpp_total -s1";...
            "%PGM%\exec\con2asc -d, -ifdnpp_sums.bin -ofdnpp_sums -s1"];


        for i = 1:length(batchCommands)
            fprintf(fid, "%s\n", batchCommands(i));
        end

        fclose(fid);
        system('Helloworld.bat');

        disp("Conc file will be read..");
        
        Conc = readtable("fdnpp_total.txt");
        Sums = readtable("fdnpp_sums.txt");
        Conc.UTC = datetime(Conc.YEAR,Conc.MO,Conc.DA,Conc.HR,0,0);
        Sums.UTC = datetime(Sums.YEAR,Sums.MO,Sums.DA,Sums.HR,0,0);
        Conc = movevars(Conc,"UTC","Before",1);
        Sums = movevars(Sums,"UTC","Before",1);
        Conc(:,["YEAR" "MO" "DA" "HR"])= [];
        Sums(:,["YEAR" "MO" "DA" "HR"])= [];

        Variable = Conc.Properties.VariableNames(4:end);
        ldata = groupsummary(Conc,"UTC",["max","mean","std"],Variable,IncludeMissingGroups=false);
        u = unique(Conc.UTC);
        [ConcCell,CentroidCell] = deal(cell(numel(u), 1));
        for d = 1:numel(u)
             ConcCell{d,:} = Conc(Conc.UTC==u(d),:);
             % RegionCell{d,:} = unique(ConcCell{d}.Region);
             % RegionCell{d,:} = nonzeros(unique(cat(2,ConcCell{d}.Region)));
             [c1,c2] = centroid(polyshape(ConcCell{d}.LAT,ConcCell{d}.LON)); %,Trajectory = Traj
             CentroidCell{d,:} = [c1,c2];
         end
        ldata.Concentration = ConcCell;
        % % ldata.RegionCell = RegionCell;
        ldata.Centroid = CentroidCell;
        % % ldata.Area = zeros(height(ldata),1);

        for i = 1:height(ldata)
            ldata.Area(i) = areaint(ldata.Concentration{i}.LAT,ldata.Concentration{i}.LON,wgs84Ellipsoid("km"));
        end

        files = dir('CG*.txt');

        files = files(~[files.isdir]);

        Conc_allData = [];

        for i = 1:length(files)
            filename = files(i).name;
            disp(['Reading ', filename]);


            data = readmatrix(filename);

            % Concatenate
            Conc_allData = [Conc_allData; data];
        end

        Conc_alldata = array2table(Conc_allData);
        Conc_alldata.UTC = datetime(Conc_alldata{:,1},Conc_alldata{:,2},Conc_alldata{:,3},Conc_alldata{:,4},0,0);
        Conc_alldata = movevars(Conc_alldata,"UTC","Before",1);
        Conc_alldata(:,2:5)=[];
        Conc_alldata.Properties.VariableNames = ["UTC","LAT","LON","Cs137000","Cs137010","I131000","I131010"];

        NewData.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(t)),month(Time(t)),day(Time(t)),hour(Time(t)))).Conc = Conc_alldata;
        NewData.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(t)),month(Time(t)),day(Time(t)),hour(Time(t)))).Total = Conc;
        NewData.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(t)),month(Time(t)),day(Time(t)),hour(Time(t)))).ldata = ldata;
        NewData.(Location(l)).(sprintf("T%02d%02d%02d%02d",year(Time(t)),month(Time(t)),day(Time(t)),hour(Time(t)))).Sums = Sums;
    end

    CompletedSimulation = CompletedSimulation+1;
    disp("Percentage completed :" + (CompletedSimulation*100)./(numel(Time)*(numel(Location)))+"%")
    Innertoc = toc(Innertic);
    disp("Time taken : ",CompletedSimulation);

end


disp('Complete.')

