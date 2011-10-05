function stereo2readable_bulk(spc, dataDir, dT)

    % 0. Constants

    Te = 13e4; % K
    k = 1.3806503e-23; % m^2*kg*s^-2*K^-1
    m = 1.67262158e-27; % kg
    
    if ~strcmp(dataDir(end), '/')
        dataDir = [dataDir, '/'];
    end
    
    dates1 = [];
    dates2 = [];
    dates3 = [];
    Np = [];
    Vp = [];
    Tp = [];
    beta = [];
    Br = [];
    Bt = [];
    Bn = [];
    B = [];
    Vr = [];
    Vt = [];
    Vn = [];
    
    yearFolders = dir(dataDir);
    
    for i = 1:max(size(yearFolders))
        
        yearFolder = yearFolders(i);
        if strcmp('.', yearFolder.name) || strcmp('..', yearFolder.name)
            continue
        end
        yearDir = [dataDir, yearFolder.name, '/'];
        
        % 1. Get plasma data
        year = str2num(yearFolder.name); %#ok<ST2NM>
        fp = fopen([yearDir, num2str(year), '_P.dat']);
        fline = 'tmp';
        hl = 0;
        while ischar(fline) && ~strcmp(fline, 'DATA:')
            fline = fgetl(fp);
            hl = hl+1;
        end
        fclose(fp);

        [years, months, days, hours, minutes, seconds, milliseconds, Np0, Vp0, Tp0, beta0] = textread([yearDir, num2str(year), '_P.dat'], '%u%u%u%u%u%u%u%f%f%f%*f%f%*f', 'headerlines', hl);
        dates10 = datenum(years, months, days, hours, minutes, seconds+milliseconds);
        
        Np = [Np; Np0]; %#ok<AGROW>
        Vp = [Vp; Vp0]; %#ok<AGROW>
        Tp = [Tp; Tp0]; %#ok<AGROW>
        beta = [beta; beta0]; %#ok<AGROW>
        dates1 = [dates1; dates10]; %#ok<AGROW>
        
        % 2. Get magnetic data
    
        fp = fopen([yearDir, num2str(year), '_M.dat']);
        fline = 'tmp';
        %spc = {};
        hl = 0;
        while ischar(fline) && ~strcmp(fline, 'DATA:')
            fline = fgetl(fp);
            %if ~ischar(spc)
            %    spc = regexp(fline, 'STEREO (A|B)', 'match');
            %    spc = cat(2, spc{:});
            %end
            hl = hl+1;
        end
        fclose(fp);
        %spc = spc(end);

        [years, months, days, hours, minutes, seconds, Br0, Bt0, Bn0, B0] = textread([yearDir, num2str(year), '_M.dat'], '%u%u%u%u%u%u%*u%f%f%f%f%*f%*f%*f%*f', 'headerlines', hl);
        dates20 = datenum(years, months, days, hours, minutes, seconds);
        
        Br = [Br; Br0]; %#ok<AGROW>
        Bt = [Bt; Bt0]; %#ok<AGROW>
        Bn = [Bn; Bn0]; %#ok<AGROW>
        B = [B; B0]; %#ok<AGROW>
        dates2 = [dates2; dates20]; %#ok<AGROW>
        
        if strcmpi(spc, 'B')
            Vr = [Vr; Vp0]; %#ok<AGROW>
            Vt = [Vt; zeros(length(dates10), 1)]; %#ok<AGROW>
            Vn = [Vn; zeros(length(dates10), 1)]; %#ok<AGROW>
            dates3 = [dates3; dates10]; %#ok<AGROW>
        else
            for month = 1:12
                files = dir([dataDir, num2str(year), '/', strcat('ST', upper(spc), '_L2_PLA_1DMax_1min_', num2str(year), num2str(month, '%02u'), '_V*.txt')]);
                if size(files, 1) > 0
                    files = sort({files.name});
                    file = files{end};
                    [years, doys, hours, minutes, seconds, milliseconds, Vr0, Vt0, Vn0] = textread([dataDir, num2str(year), '/', file], '%u%u%u%u%u%f%*s%*s%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%f%f%f%*f%*f%*d%*d%*d%*d%*d%*d%*d%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f', 'headerlines', 2);
                    monthdates = datevec(doy2date(doys, years));
                    months = monthdates(:, 2);
                    days = monthdates(:, 3);
                    dates3 = [dates3; datenum([years, months, days, hours, minutes, seconds+milliseconds])]; %#ok<AGROW>
                    Vr = [Vr; Vr0]; %#ok<AGROW>
                    Vt = [Vt; Vt0]; %#ok<AGROW>
                    Vn = [Vn; Vn0]; %#ok<AGROW>
                end
            end
        end
    end
    
    [dates1, indices] = unique(dates1);
    Np = Np(indices);
    Vp = Vp(indices);
    Tp = Tp(indices);
    beta = beta(indices);
    [dates2, indices] = unique(dates2);
    Br = Br(indices);
    Bt = Bt(indices);
    Bn = Bn(indices);
    B = B(indices);
    [dates3, indices] = unique(dates3);
    Vr = Vr(indices);
    Vt = Vt(indices);
    Vn = Vn(indices);
    
    v = Np < 1e10;
    Np = Np(v);
    dates_Np = dates1(v);
    v = Vp < 1e10;
    Vp = Vp(v);
    dates_Vp = dates1(v);
    v = Tp > 0 & Tp < 1e10;
    Tp = Tp(v);
    dates_Tp = dates1(v);
    v = beta > 0 & beta < 1e10;
    beta = beta(v);
    dates_beta = dates1(v);
    v = B < 1e10;
    B = B(v);
    dates_B= dates2(v);
    v = Br < 1e10;
    Br = Br(v);
    dates_Br= dates2(v);
    v = Bt < 1e10;
    Bt = Bt(v);
    dates_Bt= dates2(v);
    v = Bn < 1e10;
    Bn = Bn(v);
    dates_Bn= dates2(v);
    v = abs(Vr) < 1e10;
    Vr = Vr(v);
    dates_Vr= dates3(v);
    v = abs(Vt) < 1e10;
    Vt = Vt(v);
    dates_Vt= dates3(v);
    v = abs(Vn) < 1e10;
    Vn = Vn(v);
    dates_Vn= dates3(v);
    
    date_start = max([dates_Np(1), dates_Vp(1), dates_Tp(1), dates_beta(1), ...
        dates_B(1), dates_Br(1), dates_Bt(1), dates_Bn(1), ...
        dates_Vr(1), dates_Vt(1), dates_Vn(1)]);
    date_start_vector = datevec(date_start);
    if date_start_vector(6) > 0
        date_start_vector(6) = 0;
        date_start = addtodate(datenum(date_start_vector), 1, 'minute');
    end
    date_end = min([dates_Np(end), dates_Vp(end), dates_Tp(end), dates_beta(end), ...
        dates_B(end), dates_Br(end), dates_Bt(end), dates_Bn(end), ...
        dates_Vr(end), dates_Vt(end), dates_Vn(end)]);
    N = floor(etime(datevec(date_end), datevec(date_start))/dT)+1;
    dt = datenum(0, 0, 0, 0, 0, dT);
    
    newdates = (0:N-1)*dt + date_start;
    
    tsc = tscollection(newdates);
    tsc = addts(tsc, resample(timeseries(Np, dates_Np, 'name', 'Np'), newdates));
    tsc = addts(tsc, resample(timeseries(Vp, dates_Vp, 'name', 'Vp'), newdates));
    tsc = addts(tsc, resample(timeseries(Tp, dates_Tp, 'name', 'Tp'), newdates));
    tsc = addts(tsc, resample(timeseries(beta, dates_beta, 'name', 'beta'), newdates));
    tsc = addts(tsc, resample(timeseries(B, dates_B, 'name', 'B'), newdates));
    tsc = addts(tsc, resample(timeseries(Br, dates_Br, 'name', 'Br'), newdates));
    tsc = addts(tsc, resample(timeseries(Bt, dates_Bt, 'name', 'Bt'), newdates));
    tsc = addts(tsc, resample(timeseries(Bn, dates_Bn, 'name', 'Bn'), newdates));
    tsc = addts(tsc, resample(timeseries(Vr, dates_Vr, 'name', 'Vr'), newdates));
    tsc = addts(tsc, resample(timeseries(Vt, dates_Vt, 'name', 'Vt'), newdates));
    tsc = addts(tsc, resample(timeseries(Vn, dates_Vn, 'name', 'Vn'), newdates));
    
    tsc = addts(tsc, timeseries(tsc.Np.data(:)*1e6.*(tsc.Tp.data(:)+ones(N, 1)*Te)*k*1e9, newdates, 'name', 'Pth'));
    tsc = addts(tsc, timeseries((2*k*tsc.Tp.data(:)/m).^0.5/1e3, newdates, 'name', 'Vth'));
    
    % 4. Save
    
    fp = fopen(strcat('./', 'stereo_', lower(spc), '_', num2str(dT), '.dat'), 'w');
    fprintf(fp, '# STEREO-%s (RTN) %s - %s (dT = %us)\n', upper(spc), datestr(newdates(1), 'yyyy-mm-dd HH:MM'), datestr(newdates(end), 'yyyy-mm-dd HH:MM'), dT);
    fprintf(fp, '#  y   m   d   H   M   S     |B|      Br      Bt      Bn     Vp      Vr      Vt      Vn      Pth      Np         Tp    Vth      beta\n');
    fprintf(fp, '%4u  %02u  %02u  %02u  %02u  %02u  %6.2f  %6.2f  %6.2f  %6.2f  %5.1f  %6.1f  %6.1f  %6.1f  %7.5f  %6.3f  %9.1f  %5.1f  %8.5f\n', [datevec(newdates), tsc.B.data(:), tsc.Br.data(:), tsc.Bt.data(:), tsc.Bn.data(:), tsc.Vp.data(:), tsc.Vr.data(:), tsc.Vt.data(:), tsc.Vn.data(:), tsc.Pth.data(:), tsc.Np.data(:), tsc.Tp.data(:), tsc.Vth.data(:), tsc.beta.data(:)]');
    fclose(fp);
    
end
