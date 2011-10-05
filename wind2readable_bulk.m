function wind2readable_bulk(dataDir, dT)

    % 0. Constants

    Te = 13e4; % K
    k = 1.3806503e-23; % m^2*kg*s^-2*K^-1
    m = 1.67262158e-27; % kg
    
    if ~strcmp(dataDir(end), '/')
        dataDir = [dataDir, '/'];
    end
    
    dates = [];
    B = [];
    Bx = [];
    By = [];
    Bz = [];
    Vp = [];
    Vx = [];
    Vy = [];
    Vz = [];
    Vth = [];
    Np = [];
    beta = [];
    
    yearFiles = dir([dataDir, '*.dat']);
    
    for i = 1:max(size(yearFiles))

        yearFile= yearFiles(i).name;
        
        % 1. Get data from file

        [years, doys, milliseconds, B0, Bx0, By0, Bz0, Vp0, Vx0, Vy0, Vz0, Vth0, Np0, beta0] = textread([dataDir, yearFile], '%u%u%u%f%f%f%f%f%f%f%f%f%f%f', 'headerlines', 3);

        dates = [dates; datenum(doy2date(doys, years)) + datenum(milliseconds/1e3*[0, 0, 0, 0, 0, 1])]; %#ok<AGROW>
        B = [B; B0]; %#ok<AGROW>
        Bx = [Bx; Bx0]; %#ok<AGROW>
        By = [By; By0]; %#ok<AGROW>
        Bz = [Bz; Bz0]; %#ok<AGROW>
        Vp = [Vp; Vp0]; %#ok<AGROW>
        Vx = [Vx; Vx0]; %#ok<AGROW>
        Vy = [Vy; Vy0]; %#ok<AGROW>
        Vz = [Vz; Vz0]; %#ok<AGROW>
        Vth = [Vth; Vth0]; %#ok<AGROW>
        Np = [Np; Np0]; %#ok<AGROW>
        beta = [beta; beta0]; %#ok<AGROW>
        
    end
    
    [dates, indices] = unique(dates);
    tsc = tscollection(dates);
    tsc = addts(tsc, timeseries(B(indices), dates, 'name', 'B'));
    tsc = addts(tsc, timeseries(Bx(indices), dates, 'name', 'Bx'));
    tsc = addts(tsc, timeseries(By(indices), dates, 'name', 'By'));
    tsc = addts(tsc, timeseries(Bz(indices), dates, 'name', 'Bz'));
    tsc = addts(tsc, timeseries(Vp(indices), dates, 'name', 'Vp'));
    tsc = addts(tsc, timeseries(Vx(indices), dates, 'name', 'Vx'));
    tsc = addts(tsc, timeseries(Vy(indices), dates, 'name', 'Vy'));
    tsc = addts(tsc, timeseries(Vz(indices), dates, 'name', 'Vz'));
    tsc = addts(tsc, timeseries(Vth(indices), dates, 'name', 'Vth'));
    tsc = addts(tsc, timeseries(Np(indices), dates, 'name', 'Np'));
    tsc = addts(tsc, timeseries(beta(indices), dates, 'name', 'beta'));

    date_start = dates(1);
    date_start_vector = datevec(date_start);
    if date_start_vector(6) > 0
        date_start_vector(6) = 0;
        date_start = addtodate(datenum(date_start_vector), 1, 'minute');
    end
    date_end = dates(end);

    N = floor(etime(datevec(date_end), datevec(date_start))/dT)+1;
    dt = datenum(0, 0, 0, 0, 0, dT);

    newdates = (0:N-1)*dt + date_start;

    tsc = resample(tsc, newdates);

    tsc = addts(tsc, timeseries((tsc.Vth.data(:)*1e3).^2*m/2/k, newdates, 'name', 'Tp'));
    tsc = addts(tsc, timeseries(tsc.Np.data(:)*1e6*k.*(tsc.Tp.data(:)+Te*ones(N, 1))*1e9, newdates, 'name', 'Pth'));

    % 2. Save

    fp = fopen(['./', 'wind_', num2str(dT), '.dat'], 'w');
    fprintf(fp, '# WIND (GSE) %s - %s (dT = %us)\n', datestr(newdates(1), 'yyyy-mm-dd HH:MM'), datestr(newdates(end), 'yyyy-mm-dd HH:MM'), dT);
    fprintf(fp, '#  y   m   d   H   M   S     |B|      Bx      By      Bz     Vp      Vx      Vy      Vz      Pth      Np         Tp    Vth      beta\n');
    fprintf(fp, '%4u  %02u  %02u  %02u  %02u  %02u  %6.2f  %6.2f  %6.2f  %6.2f  %5.1f  %6.1f  %6.1f  %6.1f  %7.5f  %6.3f  %9.1f  %5.1f  %8.5f\n', [datevec(newdates), tsc.B.data(:), tsc.Bx.data(:), tsc.By.data(:), tsc.Bz.data(:), tsc.Vp.data(:), tsc.Vx.data(:), tsc.Vy.data(:), tsc.Vz.data(:), tsc.Pth.data(:), tsc.Np.data(:), tsc.Tp.data(:), tsc.Vth.data(:), tsc.beta.data(:)]');
    fclose(fp);
end
