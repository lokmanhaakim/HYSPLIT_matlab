%% Open ftp setting
ftpobj = ftp("ftp://arlftp.arlhq.noaa.gov");
%% Check directory 
dir(ftpobj,"archives/gdas1")
%% Change directory
cd(ftpobj,"archives/gdas1")
%% Download file
for years = [23]
    for month = ["may"]
        for week = 1
            gdasstr = sprintf("gdas1.%s%d.w%d",month,years,week); 
            mget(ftpobj,gdasstr);
        end
    end
end

%% Extract throguh looping
years = [11]; % Define years as an array if you want to iterate over multiple years
months = ["mar" "apr"]; % First 3 letter
week = 2:5; % Assuming you want to download files for week 1

for year = years
    for monthIndex = 1:length(months)
        for week = week
        month = months(monthIndex);
        gdasstr = sprintf('gdas1.%s%d.w%d', month, year, week); 
        
        % Check if the file already exists in the current directory
        fileExists = any(strcmp({dir('*').name}, gdasstr));
        
        % If the file doesn't exist, then perform mget
        if ~fileExists
            mget(ftpobj, gdasstr);
        else
            disp(['File ', gdasstr, ' already exists in the current directory.']);
        end

        end
    end
end
