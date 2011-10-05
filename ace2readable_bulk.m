function ace2readable_bulk(dataDir, dT)

    % 0. Constants

    Te = 13e4; % K
    k = 1.3806503e-23; % m^2*kg*s^-2*K^-1
    m = 1.67262158e-27; % kg
    mu = 1.25e-6; % Tm/A
    
    if ~strcmp(dataDir(end), '/')
        dataDir = [dataDir, '/'];
    end
    
    dates1 = [];
    dates2 = [];
    Np = [];
    Vp = [];
    Tp = [];
    Bx = [];
    By = [];
    Bz = [];
    B = [];
    Vx = [];
    Vy = [];
    Vz = [];
    
    yearFolders = dir(dataDir);
    
    for i = 1:max(size(yearFolders))
        
        yearFolder = yearFolders(i);
        if strcmp('.', yearFolder.name) || strcmp('..', yearFolder.name)
            continue
        end
        yearDir = [dataDir, yearFolder.name, '/'];
        
        % 1. Get plasma data
        year = str2num(yearFolder.name); %#ok<ST2NM>
        fp = fopen([yearDir, num2str(year), '_P.txt']);
        fline = 'tmp';
        hl = 0;
        while ischar(fline) && ~strcmp(fline, 'BEGIN DATA')
            fline = fgetl(fp);
            hl = hl+1;
        end
        fclose(fp);

        [years, doys, hours, minutes, seconds, Np0, Tp0, Vp0, Vx0, Vy0, Vz0] = ...
            textread([yearDir, num2str(year), '_P.txt'], '%u%u%u%u%f%f%f%f%f%f%f', 'headerlines', hl);
        [~, months, days, ~, ~, ~] = datevec(doy2date(doys, years));
        dates10 = datenum(years, months, days, hours, minutes, seconds);
        
        Np = [Np; Np0]; %#ok<AGROW>
        Tp = [Tp; Tp0]; %#ok<AGROW>
        Vp = [Vp; Vp0]; %#ok<AGROW>
        Vx = [Vx; Vx0]; %#ok<AGROW>
        Vy = [Vy; Vy0]; %#ok<AGROW>
        Vz = [Vz; Vz0]; %#ok<AGROW>
        dates1 = [dates1; dates10]; %#ok<AGROW>
        
        % 2. Get magnetic data
    
        fp = fopen([yearDir, num2str(year), '_M.txt']);
        fline = 'tmp';
        hl = 0;
        while ischar(fline) && ~strcmp(fline, 'BEGIN DATA')
            fline = fgetl(fp);
            hl = hl+1;
        end
        fclose(fp);

        [years, doys, hours, minutes, seconds, B0, Bx0, By0, Bz0] = ...
            textread([yearDir, num2str(year), '_M.txt'], '%u%u%u%u%f%f%f%f%f', 'headerlines', hl);
        [~, months, days, ~, ~, ~] = datevec(doy2date(doys, years));
        dates20 = datenum(years, months, days, hours, minutes, seconds);
        
        B = [B; B0]; %#ok<AGROW>
        Bx = [Bx; Bx0]; %#ok<AGROW>
        By = [By; By0]; %#ok<AGROW>
        Bz = [Bz; Bz0]; %#ok<AGROW>
        dates2 = [dates2; dates20]; %#ok<AGROW>
    end
    
    [dates1, indices] = unique(dates1);
    Np = Np(indices);
    Tp = Tp(indices);
    Vp = Vp(indices);
    Vx = Vx(indices);
    Vy = Vy(indices);
    Vz = Vz(indices);
    [dates2, indices] = unique(dates2);
    B = B(indices);
    Bx = Bx(indices);
    By = By(indices);
    Bz = Bz(indices);
    
    v = Np > 0 & Np < 1e10;
    Np = Np(v);
    dates_Np = dates1(v);
    v = Tp > 0 & Tp < 1e10;
    Tp = Tp(v);
    dates_Tp = dates1(v);
    v = Vp > -999 & Vp < 1e10;
    Vp = Vp(v);
    dates_Vp = dates1(v);
    v = Vx > -999 & Vx < 1e10;
    Vx = Vx(v);
    dates_Vx = dates1(v);
    v = Vy > -999 & Vy < 1e10;
    Vy = Vy(v);
    dates_Vy = dates1(v);
    v = Vz > -999 & Vz < 1e10;
    Vz = Vz(v);
    dates_Vz = dates1(v);
    v = B > -999 & B < 1e10;
    B = B(v);
    dates_B= dates2(v);
    v = Bx > -999 & Bx < 1e10;
    Bx = Bx(v);
    dates_Bx= dates2(v);
    v = By > -999 & By < 1e10;
    By = By(v);
    dates_By= dates2(v);
    v = Bz > -999 & Bz < 1e10;
    Bz = Bz(v);
    dates_Bz= dates2(v);
    
    date_start = max([dates_Np(1), dates_Tp(1), dates_Vp(1), ...
        dates_Vx(1), dates_Vy(1), dates_Vz(1), ...
        dates_B(1), dates_Bx(1), dates_By(1), dates_Bz(1)]);
    date_start_vector = datevec(date_start);
    if date_start_vector(6) > 0
        date_start_vector(6) = 0;
        date_start = addtodate(datenum(date_start_vector), 1, 'minute');
    end
    date_end = min([dates_Np(end), dates_Tp(end), dates_Vp(end), ...
        dates_Vx(end), dates_Vy(end), dates_Vz(end), ...
        dates_B(end), dates_Bx(end), dates_By(end), dates_Bz(end)]);
    N = floor(etime(datevec(date_end), datevec(date_start))/dT)+1;
    dt = datenum(0, 0, 0, 0, 0, dT);
    
    newdates = (0:N-1)*dt + date_start;
    
    tsc = tscollection(newdates);
    tsc = addts(tsc, resample(timeseries(Np, dates_Np, 'name', 'Np'), newdates));
    tsc = addts(tsc, resample(timeseries(Tp, dates_Tp, 'name', 'Tp'), newdates));
    tsc = addts(tsc, resample(timeseries(Vp, dates_Vp, 'name', 'Vp'), newdates));
    tsc = addts(tsc, resample(timeseries(Vx, dates_Vx, 'name', 'Vx'), newdates));
    tsc = addts(tsc, resample(timeseries(Vy, dates_Vy, 'name', 'Vy'), newdates));
    tsc = addts(tsc, resample(timeseries(Vz, dates_Vz, 'name', 'Vz'), newdates));
    tsc = addts(tsc, resample(timeseries(B, dates_B, 'name', 'B'), newdates));
    tsc = addts(tsc, resample(timeseries(Bx, dates_Bx, 'name', 'Bx'), newdates));
    tsc = addts(tsc, resample(timeseries(By, dates_By, 'name', 'By'), newdates));
    tsc = addts(tsc, resample(timeseries(Bz, dates_Bz, 'name', 'Bz'), newdates));
    
    tsc = addts(tsc, timeseries(tsc.Np.data(:)*1e6.*(tsc.Tp.data(:)+ones(N, 1)*Te)*k*1e9, newdates, 'name', 'Pth'));
    tsc = addts(tsc, timeseries((2*k*tsc.Tp.data(:)/m).^0.5/1e3, newdates, 'name', 'Vth'));
    tsc = addts(tsc, timeseries(tsc.Pth.data(:)./(tsc.B.data(:)/2/mu*1e-9), newdates, 'name', 'beta'));
    
    % 4. Save
    
    fp = fopen(strcat('./', 'ace_', num2str(dT), '.dat'), 'w');
    fprintf(fp, '# ACE (GSE) %s - %s (dT = %us)\n', datestr(newdates(1), 'yyyy-mm-dd HH:MM'), datestr(newdates(end), 'yyyy-mm-dd HH:MM'), dT);
    fprintf(fp, '#  y   m   d   H   M   S     |B|      Bx      By      Bz     Vp      Vx      Vy      Vz      Pth      Np         Tp    Vth      beta\n');
    fprintf(fp, '%4u  %02u  %02u  %02u  %02u  %02u  %6.2f  %6.2f  %6.2f  %6.2f  %5.1f  %6.1f  %6.1f  %6.1f  %7.5f  %6.3f  %9.1f  %5.1f  %8.5f\n', [datevec(newdates), tsc.B.data(:), tsc.Bx.data(:), tsc.By.data(:), tsc.Bz.data(:), tsc.Vp.data(:), tsc.Vx.data(:), tsc.Vy.data(:), tsc.Vz.data(:), tsc.Pth.data(:), tsc.Np.data(:), tsc.Tp.data(:), tsc.Vth.data(:), tsc.beta.data(:)]');
    fclose(fp);
    
end
