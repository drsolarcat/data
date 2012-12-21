function themis2readable_bulk(dataDir, dT)

    % 0. Constants

    k = 1.3806503e-23; % m^2*kg*s^-2*K^-1
    m = 1.67262158e-27; % kg
    mu = 1.25e-6; % Tm/A
    
    if ~strcmp(dataDir(end), '/')
        dataDir = [dataDir, '/'];
    end
    
    % 1. Get data
    
    [dates, times, B, ~, Bx, By, Bz] = ...
        textread([dataDir '/magnetic.dat'], '%s%s%f%f%f%f%f');
    for i = 1:length(dates)
        dates1(i) = datenum([dates{i}, ' ', times{i}], 'dd-mm-yyyy HH:MM:SS.FFF');
    end
    
    v = ~isnan(B);
    B = B(v);
    dates_B = dates1(v);
    v = ~isnan(Bx);
    Bx = Bx(v);
    dates_Bx = dates1(v);
    v = ~isnan(By);
    By = By(v);
    dates_By = dates1(v);
    v = ~isnan(Bz);
    Bz = Bz(v);
    dates_Bz = dates1(v);
    
    [dates, times, Np, Pp, Vx, Vy, Vz] = ...
        textread([dataDir '/ions.dat'], '%s%s%f%f%f%f%f');
    for i = 1:length(dates)
        dates2(i) = datenum([dates{i}, ' ', times{i}], 'dd-mm-yyyy HH:MM:SS.FFF');
    end
    
    v = ~isnan(Np);
    Np = Np(v);
    dates_Np = dates2(v);
    v = ~isnan(Pp);
    Pp = Pp(v);
    dates_Pp = dates2(v);
    v = ~isnan(Vx);
    Vx = Vx(v);
    dates_Vx = dates2(v);
    v = ~isnan(Vy);
    Vy = Vy(v);
    dates_Vy = dates2(v);
    v = ~isnan(Vz);
    Vz = Vz(v);
    dates_Vz = dates2(v);
    
    [dates, times, Ne, Pe, ~, ~, ~] = ...
        textread([dataDir '/electrons.dat'], '%s%s%f%f%f%f%f');
    for i = 1:length(dates)
        dates3(i) = datenum([dates{i}, ' ', times{i}], 'dd-mm-yyyy HH:MM:SS.FFF');
    end

    v = ~isnan(Ne);
    Ne = Ne(v);
    dates_Ne = dates3(v);
    v = ~isnan(Pe);
    Pe = Pe(v);
    dates_Pe = dates3(v);
    
    dateBegin = max([dates_B(1), dates_Bx(1), dates_By(1), dates_Bz(1), ...
                     dates_Np(1), dates_Pp(1), ...
                     dates_Vx(1), dates_Vy(1), dates_Vz(1), ...
                     dates_Ne(1), dates_Pe(1)]);
    dateEnd = min([dates_B(end), ...
                   dates_Bx(end), dates_By(end), dates_Bz(end), ...
                   dates_Np(end), dates_Pp(end), ...
                   dates_Vx(end), dates_Vy(end), dates_Vz(end), ...
                   dates_Ne(end), dates_Pe(end)]);
    dt = datenum(0,0,0,0,0,0.25);
    N = floor((dateEnd-dateBegin)/dt);
    dates = (0:N-1)*dt+dateBegin;
    
    tsc = tscollection(dates);
    tsc = addts(tsc, resample(timeseries(B, dates_B, 'name', 'B'), dates)); % nT
    tsc = addts(tsc, resample(timeseries(Bx, dates_Bx, 'name', 'Bx'), dates)); % nT
    tsc = addts(tsc, resample(timeseries(By, dates_By, 'name', 'By'), dates)); % nT
    tsc = addts(tsc, resample(timeseries(Bz, dates_Bz, 'name', 'Bz'), dates)); % nT
    tsc = addts(tsc, resample(timeseries(Np, dates_Np, 'name', 'Np'), dates)); % cm^-3
    tsc = addts(tsc, resample(timeseries(Pp*1.60217646e-19*1e6*1e9, dates_Pp, 'name', 'Pp'), dates)); % nPa
    tsc = addts(tsc, resample(timeseries(Vx, dates_Vx, 'name', 'Vx'), dates)); % km/s
    tsc = addts(tsc, resample(timeseries(Vy, dates_Vy, 'name', 'Vy'), dates)); % km/s
    tsc = addts(tsc, resample(timeseries(Vz, dates_Vz, 'name', 'Vz'), dates)); % km/s
    tsc = addts(tsc, resample(timeseries(Ne, dates_Ne, 'name', 'Ne'), dates)); % cm^-3
    tsc = addts(tsc, resample(timeseries(Pe*1.60217646e-19*1e6*1e9, dates_Pe, 'name', 'Pe'), dates)); % nPa
    
    tsc = addts(tsc, timeseries(sqrt(tsc.Vx.data(:).^2+ ...
                                     tsc.Vy.data(:).^2+ ...
                                     tsc.Vz.data(:).^2), dates, 'name', 'Vp')); % km/s
                                 
    tsc = addts(tsc, timeseries(tsc.Pp.data(:)+tsc.Pe.data(:), dates, 'name', 'Pth')); % nPa
    
    tsc = addts(tsc, timeseries(tsc.Pp.data(:)*1e-9./tsc.Np.data(:)/1e6/k, dates, 'name', 'Tp')); % K
    
    tsc = addts(tsc, timeseries(sqrt(2*k*tsc.Tp.data(:)/m)/1e3, dates, 'name', 'Vth')); % km/s

    tsc = addts(tsc, timeseries(tsc.Pth.data(:)./(tsc.B.data(:)/2/mu*1e-9), dates, 'name', 'beta'));
    
    % 2. Save
    
    fp = fopen(strcat('./icme/res/', 'themis_', num2str(dT), '.dat'), 'w');
    fprintf(fp, '# THEMIS (GSM) %s - %s (dT = %us)\n', datestr(dates(1), 'yyyy-mm-dd HH:MM'), datestr(dates(end), 'yyyy-mm-dd HH:MM'), dT);
    fprintf(fp, '#  y   m   d   H   M      S     |B|      Bx      By      Bz     Vp      Vx      Vy      Vz      Pth      Np         Tp    Vth      beta\n');
    fprintf(fp, '%4u  %02u  %02u  %02u  %02u  %05.2f  %6.2f  %6.2f  %6.2f  %6.2f  %5.1f  %6.1f  %6.1f  %6.1f  %7.5f  %6.3f  %9.1f  %5.1f  %8.5f\n', [datevec(dates), tsc.B.data(:), tsc.Bx.data(:), tsc.By.data(:), tsc.Bz.data(:), tsc.Vp.data(:), tsc.Vx.data(:), tsc.Vy.data(:), tsc.Vz.data(:), tsc.Pth.data(:), tsc.Np.data(:), tsc.Tp.data(:), tsc.Vth.data(:), tsc.beta.data(:)]');
    fclose(fp);
    
end
