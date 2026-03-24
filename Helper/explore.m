classdef explore
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, 
% and Piers R. F. Barnes, Imperial College London
%
    methods (Static)

        %--------------------------------------------------------------
        % 2-PARAMETER EXPLORATION
        %--------------------------------------------------------------
        function exsol = explore2par(par_base, parnames, parvalues, JVpnts)
            % EXPLORE2PAR: explore 2 parameters using the parallel
            % computing toolbox.
            %
            % PAR_BASE   base parameter struct
            % PARNAMES   cell array with 2 parameter names
            % PARVALUES  cell array with 2 numeric vectors
            % JVPnts     number of JV points

            tic
            disp('Starting 2-parameter exploration');
            disp(['Parameter 1: ', parnames{1}]);
            disp(['Parameter 2: ', parnames{2}]);

            parval1 = cell2mat(parvalues(1));
            parval2 = cell2mat(parvalues(2));
            str1    = char(parnames(1));
            str2    = char(parnames(2));

            n1 = length(parval1);
            n2 = length(parval2);

            errorlog = zeros(n1, n2);

            % Preallocate (2D in parameters, 3D for JV/PL traces)
            A  = zeros(n1, n2);              % Voc_f
            B  = zeros(n1, n2);              % Voc_r
            C  = zeros(n1, n2);              % Jsc_f
            D  = zeros(n1, n2);              % Jsc_r
            E  = zeros(n1, n2);              % mpp_f
            F  = zeros(n1, n2);              % mpp_r
            G  = zeros(n1, n2);              % FF_f
            H  = zeros(n1, n2);              % FF_r
            J  = zeros(n1, n2, JVpnts);      % Voc_stable
            K  = zeros(n1, n2, JVpnts);      % PLint

            AA = zeros(n1, n2, JVpnts);      % Vapp_f
            BB = zeros(n1, n2, JVpnts);      % J_f
            CC = zeros(n1, n2, JVpnts);      % Vapp_r
            DD = zeros(n1, n2, JVpnts);      % J_r

            EE = zeros(n1, n2);              % n_av
            FF = zeros(n1, n2);              % p_av

            parfor i = 1:n1

                par_i = par_base;
                par_i = explore.helper(par_i, str1, parval1(i));

                if strmatch('d', parnames(1)) ~= 0
                    % sets number of points in active layer to be #nm/2
                    points = round((parval1(i)*1e7)/2);
                    if points < 100
                        points = 100;
                    end
                    par_i = explore.helper(par_i, ['layer_points', str1(2:end)], points);
                end

                % Refresh device
                par_i = refresh_device(par_i);

                Voc_f     = zeros(1, n2);
                Voc_r     = zeros(1, n2);
                Jsc_f     = zeros(1, n2);
                Jsc_r     = zeros(1, n2);
                mpp_f     = zeros(1, n2);
                mpp_r     = zeros(1, n2);
                FF_f      = zeros(1, n2);
                FF_r      = zeros(1, n2);
                n_av_vec  = zeros(1, n2);
                p_av_vec  = zeros(1, n2);
                errortemp = zeros(1, n2);

                Voc_stable = zeros(n2, JVpnts);
                PLint      = zeros(n2, JVpnts);
                Vapp_f_loc = zeros(n2, JVpnts);
                J_f_loc    = zeros(n2, JVpnts);
                Vapp_r_loc = zeros(n2, JVpnts);
                J_r_loc    = zeros(n2, JVpnts);

                for j = 1:n2
                    sol0 = []; %#ok<NASGU> % predeclare for parfor analyser
                    try
                        disp(['Run no. ', num2str((i-1)*n2 + j), ...
                              ', ', str1 ,' = ', num2str(parval1(i)), ...
                              ' , ', str2,' = ', num2str(parval2(j))]);

                        par = par_i;  % local copy for this j

                        % If second parameter is Intensity use PARVAL2,
                        % else set parameter and use Int=1.
                        if strmatch('Int', parnames(2)) == 1
                            Int = parval2(j);
                        else
                            par = explore.helper(par, str2, parval2(j));
                            Int = 1;
                            par = refresh_device(par);
                        end

                        % Equilibrium
                        soleq = equilibrate(par);

                        % Choose ionic or purely electronic solution
                        if isfield(soleq, 'ion')
                            sol0 = soleq.ion;
                        elseif isfield(soleq, 'el')
                            sol0 = soleq.el;
                        else
                            error('Equilibrium solution has neither .ion nor .el field.');
                        end

                        % JV from 0–1.2 V @ 50 mV/s
                        JV    = doJV(sol0, 50e-3, JVpnts, Int, 1, 0, 1.2, 2);
                        stats = dfana.JVstats(JV);

                        % Store JV traces
                        Vapp_f_loc(j,:) = dfana.calcVapp(JV.ill.f);
                        J_f_loc(j,:)    = explore.getJtot(JV.ill.f);
                        Vapp_r_loc(j,:) = dfana.calcVapp(JV.ill.r);
                        J_r_loc(j,:)    = explore.getJtot(JV.ill.r);

                        % Scalar stats
                        Voc_f(j) = stats.Voc_f;
                        Voc_r(j) = stats.Voc_r;
                        Jsc_f(j) = stats.Jsc_f;
                        Jsc_r(j) = stats.Jsc_r;
                        mpp_f(j) = stats.mpp_f;
                        mpp_r(j) = stats.mpp_r;
                        FF_f(j)  = stats.FF_f;
                        FF_r(j)  = stats.FF_r;

                        % Steady-state Voc
                        sol_Voc         = lightonRs(sol0, 1, -10, 1, 1e6, JVpnts);
                        Voc_stable(j,:) = dfana.calcDeltaQFL(sol_Voc);
                        PLint(j,:)      = dfana.calcPLt(sol_Voc);

                        % Robust active-region averaging
                        Npts = size(sol_Voc.u, 2);
                        [i1, i2] = explore.get_active_range(par, Npts);
                        n_av_vec(j) = mean(sol_Voc.u(end, i1:i2, 1));
                        p_av_vec(j) = mean(sol_Voc.u(end, i1:i2, 2));

                        errortemp(j) = 0;

                    catch ME
                        warning(['DRIFTFUSION FAILURE: Run no. ', ...
                            num2str((i-1)*n2 + j), ', ', str1, '= ', ...
                            num2str(parval1(i)), ', ', str2, '= ', ...
                            num2str(parval2(j)), ' (', ME.message, ')']);
                        errortemp(j) = 1;
                    end
                end

                % Collect into global arrays
                A(i,:)   = Voc_f;
                B(i,:)   = Voc_r;
                C(i,:)   = Jsc_f;
                D(i,:)   = Jsc_r;
                E(i,:)   = mpp_f;
                F(i,:)   = mpp_r;
                G(i,:)   = FF_f;
                H(i,:)   = FF_r;
                J(i,:,:) = Voc_stable;
                K(i,:,:) = PLint;
                AA(i,:,:) = Vapp_f_loc;
                BB(i,:,:) = J_f_loc;
                CC(i,:,:) = Vapp_r_loc;
                DD(i,:,:) = J_r_loc;
                EE(i,:)   = n_av_vec;
                FF(i,:)   = p_av_vec;
                errorlog(i,:) = errortemp;
            end

            % Pack into output struct
            exsol.stats.Voc_f      = A;
            exsol.stats.Voc_r      = B;
            exsol.stats.Jsc_f      = C;
            exsol.stats.Jsc_r      = D;
            exsol.stats.mpp_f      = E;
            exsol.stats.mpp_r      = F;
            exsol.stats.FF_f       = G;
            exsol.stats.FF_r       = H;
            exsol.Vapp_f           = AA;
            exsol.J_f              = BB;
            exsol.Vapp_r           = CC;
            exsol.J_r              = DD;
            exsol.stats.Voc_stable = J;
            exsol.stats.PLint      = K;
            exsol.parnames         = parnames;
            exsol.parvalues        = parvalues;
            exsol.parval1          = parval1;
            exsol.parval2          = parval2;
            exsol.par_base         = par_base;
            exsol.stats.n_av       = EE;
            exsol.stats.p_av       = FF;
            exsol.errorlog         = errorlog;

            toc
        end

        %--------------------------------------------------------------
        % 1-PARAMETER EXPLORATION
        %--------------------------------------------------------------
        function exsol = explore1par(par_base, parname, parvalues, JVpnts)
        % EXPLORE1PAR: explore 1 parameter using the parallel computing
        % toolbox, but store results in 2-parameter shapes with a dummy
        % second dimension so all plot functions work.

            tic
            disp('Starting 1-parameter exploration');

            if iscell(parname)
                parname = parname{1};
            end
            str1 = char(parname);
            disp(['Parameter: ', str1]);

            if iscell(parvalues)
                parval = cell2mat(parvalues(1));
            else
                parval = parvalues;
            end

            nvals = length(parval);
            n2    = 1; % dummy second parameter dimension

            % Preallocate like 2-parameter case with n2 = 1
            A  = zeros(nvals, n2);           % Voc_f
            B  = zeros(nvals, n2);           % Voc_r
            C  = zeros(nvals, n2);           % Jsc_f
            D  = zeros(nvals, n2);           % Jsc_r
            E  = zeros(nvals, n2);           % mpp_f
            F  = zeros(nvals, n2);           % mpp_r
            G  = zeros(nvals, n2);           % FF_f
            H  = zeros(nvals, n2);           % FF_r
            Jv = zeros(nvals, n2, JVpnts);   % Voc_stable
            Kv = zeros(nvals, n2, JVpnts);   % PLint

            AA = zeros(nvals, n2, JVpnts);   % Vapp_f
            BB = zeros(nvals, n2, JVpnts);   % J_f
            CC = zeros(nvals, n2, JVpnts);   % Vapp_r
            DD = zeros(nvals, n2, JVpnts);   % J_r

            EE  = zeros(nvals, n2);          % n_av
            FFv = zeros(nvals, n2);          % p_av

            errorlog = zeros(nvals, n2);

            parfor i = 1:nvals

                par = par_base;

                % Set parameter value unless it's Int*
                if strmatch('Int', str1) == 1
                    % leave 'par' unchanged; Int set below
                else
                    par = explore.helper(par, str1, parval(i));
                end

                % If thickness parameter starting with 'd', adjust layer_points
                if strmatch('d', str1) ~= 0
                    points = round((parval(i)*1e7)/2);
                    if points < 100
                        points = 100;
                    end
                    par = explore.helper(par, ['layer_points', str1(2:end)], points);
                end

                par = refresh_device(par);

                Voc_f = 0;
                Voc_r = 0;
                Jsc_f = 0;
                Jsc_r = 0;
                mpp_f = 0;
                mpp_r = 0;
                FF_f  = 0;
                FF_r  = 0;
                n_av  = 0;
                p_av  = 0;

                Voc_stable = zeros(1, JVpnts);
                PLint      = zeros(1, JVpnts);
                Vapp_f     = zeros(1, JVpnts);
                J_f        = zeros(1, JVpnts);
                Vapp_r     = zeros(1, JVpnts);
                J_r        = zeros(1, JVpnts);

                sol0 = []; %#ok<NASGU>
                try
                    disp(['Run no. ', num2str(i), ', ', str1, ' = ', num2str(parval(i))]);

                    if strmatch('Int', str1) == 1
                        Int = parval(i);
                    else
                        Int = 1;
                    end

                    % Equilibrium
                    soleq = equilibrate(par);

                    if isfield(soleq, 'ion')
                        sol0 = soleq.ion;
                    elseif isfield(soleq, 'el')
                        sol0 = soleq.el;
                    else
                        error('Equilibrium solution has neither .ion nor .el field.');
                    end

                    % JV from 0–1.2 V @ 50 mV/s
                    JV    = doJV(sol0, 50e-3, JVpnts, Int, 1, 0, 1.2, 2);
                    stats = dfana.JVstats(JV);

                    % JV traces
                    Vapp_f(:) = dfana.calcVapp(JV.ill.f);
                    J_f(:)    = explore.getJtot(JV.ill.f);
                    Vapp_r(:) = dfana.calcVapp(JV.ill.r);
                    J_r(:)    = explore.getJtot(JV.ill.r);

                    % Scalar stats
                    Voc_f = stats.Voc_f;
                    Voc_r = stats.Voc_r;
                    Jsc_f = stats.Jsc_f;
                    Jsc_r = stats.Jsc_r;
                    mpp_f = stats.mpp_f;
                    mpp_r = stats.mpp_r;
                    FF_f  = stats.FF_f;
                    FF_r  = stats.FF_r;

                    % Steady-state Voc
                    sol_Voc       = lightonRs(sol0, 1, -10, 1, 1e6, JVpnts); 
                    Voc_stable(:) = dfana.calcDeltaQFL(sol_Voc);
                    PLint(:)      = dfana.calcPLt(sol_Voc);

                    % Robust active-region averaging
                    Npts = size(sol_Voc.u, 2);
                    [i1, i2] = explore.get_active_range(par, Npts);
                    n_av     = mean(sol_Voc.u(end, i1:i2, 1));
                    p_av     = mean(sol_Voc.u(end, i1:i2, 2));

                    errflag = 0;

                catch ME
                    warning(['DRIFTFUSION FAILURE: Run no. ', num2str(i), ...
                        ', ', str1, ' = ', num2str(parval(i)), ...
                        ' (', ME.message, ')']);
                    errflag = 1;
                end

                % Store for this i (2-param shape with n2 = 1)
                A(i,1)    = Voc_f;
                B(i,1)    = Voc_r;
                C(i,1)    = Jsc_f;
                D(i,1)    = Jsc_r;
                E(i,1)    = mpp_f;
                F(i,1)    = mpp_r;
                G(i,1)    = FF_f;
                H(i,1)    = FF_r;
                Jv(i,1,:) = Voc_stable;
                Kv(i,1,:) = PLint;
                AA(i,1,:) = Vapp_f;
                BB(i,1,:) = J_f;
                CC(i,1,:) = Vapp_r;
                DD(i,1,:) = J_r;
                EE(i,1)   = n_av;
                FFv(i,1)  = p_av;
                errorlog(i,1) = errflag;

            end % parfor

            % Pack into output struct (2-param compatible)
            exsol.stats.Voc_f       = A;
            exsol.stats.Voc_r       = B;
            exsol.stats.Jsc_f       = C;
            exsol.stats.Jsc_r       = D;
            exsol.stats.mpp_f       = E;
            exsol.stats.mpp_r       = F;
            exsol.stats.FF_f        = G;
            exsol.stats.FF_r        = H;
            exsol.Vapp_f            = AA;
            exsol.J_f               = BB;
            exsol.Vapp_r            = CC;
            exsol.J_r               = DD;
            exsol.stats.Voc_stable  = Jv;
            exsol.stats.PLint       = Kv;

            exsol.par_base          = par_base;
            exsol.stats.n_av        = EE;
            exsol.stats.p_av        = FFv;
            exsol.errorlog          = errorlog;

            % 2-parameter style metadata (dummy second param)
            exsol.parval1           = parval;
            exsol.parval2           = 1;               % dummy
            exsol.parnames          = {str1, 'dummy'};
            exsol.parvalues         = {parval, 1};

            toc
        end

        %--------------------------------------------------------------
        % Helper: total current
        %--------------------------------------------------------------
        function Jtot = getJtot(sol)
           % Extract the total current for current structure
           J = dfana.calcJ(sol, "sub");
           Jtot = J.tot(:,end);
        end

        %--------------------------------------------------------------
        % Helper: writevar (unused right now but kept)
        %--------------------------------------------------------------
        function var = writevar(var, i, j, xx, arr)
            % Allows variable-length vectors to be stored in a 3D array
            var(i, j, 1:length(xx)) = arr;
        end

        %--------------------------------------------------------------
        % Helper: set parameter in par
        %--------------------------------------------------------------
        function par = helper(par, parname, parvalue)
            eval(['par.',parname,' = parvalue;']);
        end

        %--------------------------------------------------------------
        % Plot functions
        %--------------------------------------------------------------

        function plotPL(exsol)
            % PL surface:
            % - 2-param: surf over (parval2, parval1)
            % - 1-param: falls back to line plot vs parval1.
            PL = exsol.stats.PLint;   % [n1 x n2 x JVpnts]
            PL_end = squeeze(PL(:,:,end));  % [n1 x n2] or [n1 x 1]

            n1 = numel(exsol.parval1);
            n2 = numel(exsol.parval2);

            figure(100); clf;
            if n1 > 1 && n2 > 1
                % True 2D surface
                surf(exsol.parval2, exsol.parval1, PL_end);
                ylabel(exsol.parnames{1})
                xlabel(exsol.parnames{2})
                set(gca,'YScale','log');
                zlabel('PL intensity [cm^{-2}s^{-1}]')
                shading interp
                colorbar
            else
                % 1D case: simple line vs parval1
                if n1 > 1
                    y = PL_end(:,1);
                    semilogy(exsol.parval1, y, '-o');
                    xlabel(exsol.parnames{1});
                else
                    y = PL_end(1,:);
                    semilogy(exsol.parval2, y, '-o');
                    xlabel(exsol.parnames{2});
                end
                ylabel('PL intensity [cm^{-2}s^{-1}]');
                grid on;
            end
        end

        function plotsurf(exsol, yproperty, xlogon, ylogon, zlogon)
            % Generic surface plot of a scalar stat vs parameter1, parameter2
            try
                eval(['y = exsol.stats.', yproperty, ';'])
            catch
                error('YPROPERTY is not a property contained in exsol.STATS');
            end

            y2 = squeeze(y);  % collapse singleton dims
            n1 = numel(exsol.parval1);
            n2 = numel(exsol.parval2);

            figure(101); clf;

            if n1 > 1 && n2 > 1
                % True 2D surface
                surf(exsol.parval2, exsol.parval1, y2)
                s1 = gca;
                ylabel(exsol.parnames{1})
                xlabel(exsol.parnames{2})
                zlabel(yproperty)
                xlim([exsol.parval2(1), exsol.parval2(end)])
                ylim([exsol.parval1(1), exsol.parval1(end)])

                if xlogon
                    set(s1,'XScale','log');
                else
                    set(s1,'XScale','linear');
                end
                if ylogon
                    set(s1,'YScale','log');
                else
                    set(s1,'YScale','linear');
                end
                shading interp
                cb = colorbar();
                if zlogon
                    cb.Ruler.Scale = 'log';
                    cb.Ruler.MinorTick = 'on';
                end
            else
                % 1D case – just plot vs the swept parameter
                hold on;
                if n1 > 1
                    x = exsol.parval1(:);
                    z = y2(:,1);
                    if xlogon && ~ylogon
                        semilogx(x, z, '-o');
                    elseif ~xlogon && ylogon
                        semilogy(x, z, '-o');
                    elseif xlogon && ylogon
                        loglog(x, z, '-o');
                    else
                        plot(x, z, '-o');
                    end
                    xlabel(exsol.parnames{1});
                else
                    x = exsol.parval2(:);
                    z = y2(1,:);
                    if xlogon && ~ylogon
                        semilogx(x, z, '-o');
                    elseif ~xlogon && ylogon
                        semilogy(x, z, '-o');
                    elseif xlogon && ylogon
                        loglog(x, z, '-o');
                    else
                        plot(x, z, '-o');
                    end
                    xlabel(exsol.parnames{2});
                end
                ylabel(yproperty);
                grid on;
            end
        end



        function plotstat_2D_parval1(exsol, yproperty, logx, logy)
            % For each parval1, plot y vs parval2
            eval(['y = exsol.stats.', yproperty, ';']);
            figure(103); clf; hold on;
            for i = 1:length(exsol.parval1)
                if logx == 0 && logy == 0
                    plot(exsol.parval2, y(i,:));
                elseif logx == 1 && logy == 0
                    semilogx(exsol.parval2, y(i,:));
                elseif logx == 0 && logy == 1
                    semilogy(exsol.parval2, y(i,:));
                elseif logx == 1 && logy == 1
                    loglog(exsol.parval2, y(i,:));
                end
            end
            xlabel(exsol.parnames{2});
            ylabel(yproperty);
            legstr = (num2str((exsol.parval1' + 60e-7) * 1e7));
            legend(legstr);
            % Only apply x-limits if there are at least 2 distinct parval2 points
            if numel(exsol.parval2) > 1 && (max(exsol.parval2) > min(exsol.parval2))
                xlim([min(exsol.parval2), max(exsol.parval2)]);
            end
            hold off;
        end

        function plotstat_2D_parval2(exsol, yproperty, logx, logy)
            figure(104); clf; hold on;

            eval(['y = exsol.stats.', yproperty,';']);

            if strmatch('Voc_stable', yproperty) == 1
                y = y(:,:,end);
            end

            if strmatch('Jsc_r', yproperty) == 1
                y = y;
            end

            x = (exsol.parval1' + 60e-7) * 1e7;

            for i = 1:length(exsol.parval2)
                if logx == 0 && logy == 0
                    plot(x, y(:,i));
                elseif logx == 1 && logy == 0
                    semilogx(x, y(:,i));
                elseif logx == 0 && logy == 1
                    semilogy(x, y(:,i));
                elseif logx == 1 && logy == 1
                    loglog(x, y(:,i));
                end
            end
            xlabel('Active layer thickness [nm]');
            ylabel(yproperty);
            legstr = (num2str(exsol.parval2'));
            legend(legstr);
            % Only apply x-limits if we have a proper range in parval1
            if numel(x) > 1 && (max(x) > min(x))
                xlim([min(x), max(x)]);
            end
            hold off;
        end

        % NOTE: The following profile/recombination/CE plots rely on
        % exsol.n_f, exsol.p_f, exsol.x etc., which are only saved if you
        % re-enable that part of explore2par / explore1par. They are left
        % essentially as in the original code.

        function plotfinalELx(exsol)
            figure(105); clf;
            for i=1:length(exsol.parval1)
                for j = 1:length(exsol.parval2)
                    sol.u(1,:,1) = exsol.Vf(i, j, :);
                    sol.u(1,:,2) = exsol.nf(i, j, :);
                    sol.u(1,:,3) = exsol.pf(i, j, :);
                    sol.u(1,:,4) = exsol.cf(i, j, :);

                    sol.t   = 0;
                    sol.x   = exsol.x(i,j,:);
                    sol.par = exsol.par_base;

                    % call your EL plotter here if desired
                end
            end
        end

        function plotprof_2D(exsol, yproperty, par1logical, par2logical, logx,logy)
            par = exsol.par_base;

            eval(['y = exsol.', yproperty,';']);
            if length(par1logical) > length(exsol.parval1)
                par1logical = par1logical(1:length(exsol.parval1));
            end

            if length(par2logical) > length(exsol.parval2)
                par2logical = par2logical(1:length(exsol.parval2));
            end

            y(y==0)        = NaN;
            exsol.x(exsol.x==0) = NaN;

            if strmatch(yproperty,'a_f')
                y = y-exsol.par_base.Nani(1);
            end

            parval1 = cell2mat(exsol.parvalues(1));
            str1 = char(exsol.parnames(1));

            figure(106); clf; hold on
            for i=1:length(exsol.parval1)
                if par1logical(i) == 1

                    if strmatch('dcell', str1) ~= 0
                        layerpoints = round(parval1(i)*1e7);
                        par = explore.helper(par, ['p', str1(2:end)], layerpoints);
                    end

                    par = refresh_device(par);

                    for j = 1:length(exsol.parval2)
                        if par2logical(j) == 1
                            xplot = squeeze(exsol.x(i, j, :).*1e7);
                            yplot = squeeze(y(i, j, :));
                            if logx && ~logy
                                semilogx(xplot, yplot);
                            elseif ~logx && logy
                                semilogy(xplot, yplot);
                            elseif logx && logy
                                loglog(xplot, yplot);
                            else
                                plot(xplot, yplot);
                            end
                        end
                    end
                end
            end
            hold off
            xlabel('Position [nm]')
            ylabel(yproperty)
        end

        function plotU(exsol, par1logical, par2logical,logx,logy)
            if length(par1logical) > length(exsol.parval1)
                par1logical = par1logical(1:length(exsol.parval1));
            end

            if length(par2logical) > length(exsol.parval2)
                par2logical = par2logical(1:length(exsol.parval2));
            end

            rHTL    = zeros(length(par1logical),length(par2logical));
            rHTLint = zeros(length(par1logical),length(par2logical));
            rbulk   = zeros(length(par1logical),length(par2logical));
            rETLint = zeros(length(par1logical),length(par2logical));
            rETL    = zeros(length(par1logical),length(par2logical));
            rtotal  = zeros(length(par1logical),length(par2logical));

            figure(107); clf; hold on
            par = exsol.par_base;

            for i=1:length(exsol.parval1)
                if par1logical(i) == 1
                    for j = 1:length(exsol.parval2)
                        if par2logical(j) == 1
                            par = explore.helper(par, exsol.parnames{1,1}, exsol.parval1(i));

                            if strmatch('dcell(1,4)', exsol.parnames(1)) ~= 0
                                pcontact = round(exsol.parval1(i)*1e7);
                                par.layer_points(1,4) = pcontact*1;
                            end

                            par = refresh_device(par);

                            dev = par.dev;

                            n = squeeze(exsol.n_f(i,j,:))';
                            p = squeeze(exsol.p_f(i,j,:))';
                            n = n(1:length(dev.B));
                            p = p(1:length(dev.B));

                            rbtb = dev.B.*(n.*p - dev.ni.^2);
                            rsrh = ((n.*p - dev.ni.^2)./((dev.taun.*(p+dev.pt)) + (dev.taup.*(n+dev.nt))));

                            r = rbtb + rsrh;

                            xplot = squeeze(exsol.x(i, j, :).*1e7);
                            xplot = xplot(1:length(dev.B));

                            rHTL(i,j)    = trapz(xplot(par.pcum0(1):par.pcum0(2)), r(par.pcum0(1):par.pcum0(2)));
                            rHTLint(i,j) = trapz(xplot(par.pcum(1):par.pcum(2)), r(par.pcum(1):par.pcum(2)));
                            rbulk(i,j)   = trapz(xplot(par.pcum(2):par.pcum(5)), r(par.pcum(2):par.pcum(5)));
                            rETLint(i,j) = trapz(xplot(par.pcum(5):par.pcum(6)), r(par.pcum(5):par.pcum(6)));
                            rETL(i,j)    = trapz(xplot(par.pcum(6):par.pcum(7)), r(par.pcum(6):par.pcum(7)));
                            rtotal(i,j)  = trapz(xplot, r);

                            if logx && ~logy
                                semilogx(xplot, r);
                            elseif ~logx && logy
                                semilogy(xplot, r);
                            elseif logx && logy
                                loglog(xplot, r);
                            else
                                plot(xplot, r);
                            end
                        end
                    end
                end
            end

            hold off
            xlabel('Position, x [nm]')
            ylabel('Recombination rate, r [cm^{-3}s^{-1}]')

            figure(114); clf; hold on
            for j=1:length(exsol.parval2)
                if par2logical(j) == 1
                    plot(exsol.parval1*1e7, rHTL(:,j),...
                         exsol.parval1*1e7, rHTLint(:,j),...
                         exsol.parval1*1e7, rbulk(:,j),...
                         exsol.parval1*1e7, rETLint(:,j),...
                         exsol.parval1*1e7, rETL(:,j),...
                         exsol.parval1*1e7, rtotal(:,j));
                end
            end
            xlabel('Active layer thickness [nm]')
            ylabel('Recombination rate [cm^{-3}s^{-1}]')
            legend('HTL', 'HTL interface', 'Bulk perovskite', ...
                   'ETL interface', 'ETL','Total device')
            hold off
        end

        function plotCE(exsol_Voc, exsol_eq, xlogon, ylogon, zlogon, normalise)
            % Charge extraction plots; unchanged in shape logic.

            par     = exsol_Voc.par_base;
            parval1 = cell2mat(exsol_Voc.parvalues(1));
            parval2 = cell2mat(exsol_Voc.parvalues(2));
            str1    = char(exsol_Voc.parnames(1));

            n_CE = zeros(length(exsol_Voc.parval1), length(exsol_Voc.parval2));
            p_CE = zeros(length(exsol_Voc.parval1), length(exsol_Voc.parval2));

            for i=1:length(exsol_Voc.parval1)
                par = explore.helper(par, exsol_Voc.parnames{1,1}, exsol_Voc.parval1(i));

                if strmatch('dcell', str1) ~= 0
                    layerpoints = round(parval1(i)*1e7);
                    par = explore.helper(par, ['p', str1(2:end)], layerpoints);
                end

                par = refresh_device(par);

                n_CEx = zeros(1,length(par.xx));
                p_CEx = zeros(1,length(par.xx));

                for j = 1:length(exsol_Voc.parval2)
                    n_Voc = squeeze(exsol_Voc.n_f(i,j,1:length(par.xx)))';
                    p_Voc = squeeze(exsol_Voc.p_f(i,j,1:length(par.xx)))';
                    n_eq  = squeeze(exsol_eq.n_f(i,1:length(par.xx)));
                    p_eq  = squeeze(exsol_eq.p_f(i,1:length(par.xx)));

                    n_CEx = n_Voc - n_eq;
                    p_CEx = p_Voc - p_eq;

                    n_CE(i,j) = trapz(par.xx, n_CEx);
                    p_CE(i,j) = trapz(par.xx, p_CEx);

                    if normalise
                        % implement if needed
                    end
                end
            end

            figure(109); clf; hold on
            for i = 1:length(exsol_Voc.parval1)
                loglog(exsol_Voc.parval2, n_CE(i,:))
            end
            hold off
            xlabel('Light intensity [Suns]')
            ylabel('integrated electron density [cm^{-2}]')

            figure(110); clf;
            d_active = round(exsol_Voc.parval1*1e7)+60;
            surf(exsol_Voc.parval2, d_active , n_CE)
            s1 = gca;
            xlabel('Light intensity [Sun]')
            ylabel('Active layer thickness [nm]')
            zlabel('integrated electron density [cm^{-2}]')
            xlim([exsol_Voc.parval2(1), exsol_Voc.parval2(end)])
            ylim([d_active(1),d_active(end)])

            if xlogon
                set(s1,'XScale','log');
            else
                set(s1,'XScale','linear');
            end
            if ylogon
                set(s1,'YScale','log');
            else
                set(s1,'YScale','linear');
            end
            shading interp
            cb = colorbar();
            if zlogon
                cb.Ruler.Scale = 'log';
                cb.Ruler.MinorTick = 'on';
            end

            figure(111); clf;
            d_active = round(exsol_Voc.parval1*1e7)+60;
            surf(exsol_Voc.parval2, d_active, p_CE)
            s1 = gca;
            xlabel('Light intensity [Sun]')
            ylabel('Active layer thickness [nm]')
            zlabel('integrated hole density [cm^{-2}]')
            xlim([exsol_Voc.parval2(1), exsol_Voc.parval2(end)])
            ylim([d_active(1),d_active(end)])

            if xlogon
                set(s1,'XScale','log');
            else
                set(s1,'XScale','linear');
            end
            if ylogon
                set(s1,'YScale','log');
            else
                set(s1,'YScale','linear');
            end
            shading interp
            cb = colorbar();
            if zlogon
                cb.Ruler.Scale = 'log';
                cb.Ruler.MinorTick = 'on';
            end
        end

        function plotJV(exsol, par1logical, par2logical)
            % Works for 1-parameter (n2=1) and 2-parameter cases
            figure(112); clf;

            n1 = length(exsol.parval1);
            n2 = length(exsol.parval2);

            if length(par1logical) > n1
                par1logical = par1logical(1:n1);
            end

            if length(par2logical) > n2
                par2logical = par2logical(1:n2);
            end

            hold on
            for i=1:n1
                if par1logical(i) == 1
                    for j = 1:n2
                        if par2logical(j) == 1
                            plot(squeeze(exsol.Vapp_f(i,j,:)), squeeze(exsol.J_f(i,j,:)))
                        end
                    end
                end
            end
            xlabel('Applied Voltage [V]')
            ylabel('Current density [A cm^{-2}]')
            ylim([-30e-3, 10e-3])
            hold off
        end

        function [i1, i2] = get_active_range(par, N)
            % Returns a safe index range [i1, i2] over which to average n and p.
            % N is the total number of spatial grid points (size(sol.u,2)).

            if isfield(par, 'pcum') && ~isempty(par.pcum)
                pc = par.pcum(:).';  % row vector
                if numel(pc) >= 5
                    % Original behaviour (e.g. 7-region perovskite stack)
                    i1 = pc(2);
                    i2 = pc(5);
                elseif numel(pc) >= 2
                    % Fall back to "everything between first and last"
                    i1 = pc(1);
                    i2 = pc(end);
                else
                    % Only one cumulative index? Just use whole device
                    i1 = 1;
                    i2 = N;
                end
            else
                % No pcum at all – average over entire device
                i1 = 1;
                i2 = N;
            end

            % Clamp to valid range and ensure i1 <= i2
            i1 = max(1, min(i1, N));
            i2 = max(1, min(i2, N));
            if i2 < i1
                tmp = i1; i1 = i2; i2 = tmp; 
            end
        end

    end
end
