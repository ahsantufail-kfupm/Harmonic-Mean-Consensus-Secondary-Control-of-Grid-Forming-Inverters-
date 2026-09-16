%% ========================================================================
%  comparison_correct.m   (MATLAB)
%  Quantitative comparison — same problem, correct implementations
%
%  All four methods solve EXACTLY the same problem:
%    distributed secondary control for frequency restoration +
%    proportional active power sharing in islanded AC microgrids.
%
%  Methods:
%   M1 – Proposed       : harmonic-mean weights, consensus on r_i = P_i/Pmax_i
%   M2 – Simpson-Porco 2015 [simpson2] : DAPI controller, Metropolis weights
%         Control law (their eq. 6b):
%         k_i * dOmega_i/dt = -(omega_i - omega*) - sum_j a_ij*(Omega_i - Omega_j)
%         setpoint: u_i = -k_i*dw_i + Omega_i   (secondary variable Omega)
%   M3 – Ding et al. 2019 [ding-event] : per-unit consensus, unweighted graph
%         Control law (their eq. 3-4):
%         u_dot_i = c1*(omega* - omega_i) + c2*sum_j a_ij*(r_j - r_i)
%         same as proposed but a_ij = 1 (unweighted Laplacian)
%   M4 – Shi et al. 2020 [shi-piconsensus] : PI-consensus on per-unit loading
%         Control law (their eq. 7-8):
%         u_dot_i = kP*(omega* - omega_i) + kI*xi_i + c*sum_j a_ij*(r_j - r_i)
%         xi_dot_i = omega* - omega_i  (integral state)
%         Metropolis weights
%
%  Test: Scenario 1 — +30 MW load step at t=2s on six-unit test system.
%  All methods: same plant, droop gains, comm. graph, update period Ts=0.02s.
%
%  Metrics:
%   - Sharing settling time (time for max(r_i)-min(r_i) < 1%)
%   - Number of secondary updates to settle
%   - Frequency nadir and peak
%   - Initial RoCoF
%   - Final sharing error at t=40s
%   - CPU time per secondary update
%
%  Outputs: fig_comparison_correct.png, tab_comparison_correct.tex
%  ========================================================================

function comparison_correct()
    clc; close all;
    fprintf('=====================================================\n');
    fprintf(' Comparative simulation: same problem, 4 methods\n');
    fprintf('=====================================================\n\n');

    event = {2.0, 'load', [4, 0.30]};   % Scenario 1 disturbance

    fprintf('Running M1: Proposed (harmonic-mean)...\n');
    R1 = run_method(event, 'proposed');
    fprintf('Running M2: Simpson-Porco 2015 (DAPI)...\n');
    R2 = run_method(event, 'simpsonporco');
    fprintf('Running M3: Ding et al. 2019 (event-triggered, unweighted)...\n');
    R3 = run_method(event, 'ding');
    fprintf('Running M4: Shi et al. 2020 (PI-consensus)...\n');
    R4 = run_method(event, 'shi');
    fprintf('\nAll done.\n\n');

    labels = {'Proposed (harmonic-mean)', ...
              'Simpson-Porco et al. [16]', ...
              'Ding et al. [20]', ...
              'Shi et al. [18]'};
    cites  = {'---','\cite{simpson2}','\cite{ding-event}','\cite{shi-piconsensus}'};
    RR = {R1, R2, R3, R4};

    %% Metrics
    p = params();
    t_step = 2.0;
    gen_idx = find(p.gen0);

    nadir   = zeros(1,4); peak    = zeros(1,4);
    rocof   = zeros(1,4); settle  = zeros(1,4);
    sharerr = zeros(1,4); niter   = zeros(1,4);
    lam2v   = zeros(1,4); cputime_us = zeros(1,4);

    for m = 1:4
        Rm = RR{m};
        tv = Rm.t;
        fCOI = mean(Rm.F(:, gen_idx), 2);
        R_act = Rm.R(:, gen_idx);
        R_act(isnan(R_act)) = 0;
        share_err = max(R_act,[],2) - min(R_act,[],2);

        idx_post = tv >= t_step;
        f_post   = fCOI(idx_post);
        t_post   = tv(idx_post);
        se_post  = share_err(idx_post);

        nadir(m) = min(f_post);
        peak(m)  = max(f_post);

        % RoCoF: average over first 200ms after step
        i0 = find(tv >= t_step, 1);
        i1 = find(tv >= t_step + 0.2, 1);
        rocof(m) = (fCOI(i1) - fCOI(i0)) / (tv(i1) - tv(i0));

        % Settling time: first sample AFTER the post-step error peak where the
        % sharing error falls below 1% (searching from t_step directly returns
        % the step instant, before the disturbance-induced error has grown)
        [~, ipk] = max(se_post);
        idx_s = find(se_post(ipk:end) < 0.01, 1);
        if ~isempty(idx_s)
            idx_s = ipk + idx_s - 1;
            settle(m) = t_post(idx_s);
            % Number of secondary updates from step to settling
            niter(m) = round((settle(m) - t_step) / p.TS);
        else
            settle(m) = NaN;
            niter(m)  = NaN;
        end

        sharerr(m)    = share_err(end);
        lam2v(m)      = Rm.lambda2;
        cputime_us(m) = Rm.cpu_us;  % microseconds per secondary update

        fprintf('%s:\n', labels{m});
        fprintf('  lambda2=%.4f  nadir=%.4f Hz  RoCoF=%.4f Hz/s\n', ...
                lam2v(m), nadir(m), rocof(m));
        fprintf('  settle=%.3f s  updates=%d  share_err=%.2e  CPU=%.1f us\n\n', ...
                settle(m), niter(m), sharerr(m), cputime_us(m));
    end

    %% Figure
    make_figure(RR, labels, lam2v, settle, niter, p, t_step);

    %% LaTeX table
    write_table(labels, cites, lam2v, nadir, rocof, settle, niter, ...
                sharerr, cputime_us);

    fprintf('Figure: fig_comparison_correct.png\n');
    fprintf('Table:  tab_comparison_correct.tex\n');
end


%% ========================================================================
function R = run_method(event_in, method)
    p  = params();
    N  = p.N;
    dt = 1e-3;   T = 40;
    steps = round(T/dt);
    nsub  = round(p.TS/dt);   % secondary period = 20 ms = 20 steps

    gen   = p.gen0;
    lines = p.LINES;
    bus   = true(1,N);
    Pload = p.Pload0;
    Pmax  = p.Pmax;  Pmin = p.Pmin;

    % Build weight matrices
    A_harm  = Amat(p.COMM, gen, 'harmonic',  N);
    A_metro = Amat(p.COMM, gen, 'metropolis', N);
    A_unw   = Amat_ones(p.COMM, gen, N);   % unweighted: a_ij=1

    W_harm  = Wperron(A_harm,  gen);
    W_metro = Wperron(A_metro, gen);
    W_unw   = Wperron(A_unw,   gen);

    switch method
        case 'proposed',      A_use=A_harm;  W_use=W_harm;
        case 'simpsonporco',  A_use=A_metro; W_use=W_metro;
        case 'ding',          A_use=A_unw;   W_use=W_unw;
        case 'shi',           A_use=A_metro; W_use=W_metro;
    end
    l2 = lam2_val(A_use, gen);

    % Initial conditions (same for all methods)
    r0 = sum(p.Pload0(gen)) / sum(Pmax(gen));
    u  = r0 * Pmax .* gen;
    delta = zeros(1,N);
    dw    = zeros(1,N);

    % Extra states
    Omega = zeros(1,N);   % Simpson-Porco DAPI integral state
    Omega(gen) = u(gen) ./ p.KDR(gen);   % consistent init: u_i = KDR_i * Omega_i
    xi    = zeros(1,N);   % Shi PI integral of freq error

    inj = @(uu,dwl) max(min(uu - p.KDR.*dwl, Pmax), Pmin) .* gen;

    % Storage
    Fmat = zeros(steps, N);
    Rmat = nan(steps, N);
    tvec = zeros(steps, 1);
    ev = 1;
    % wrap event for indexing
    if ~iscell(event_in{1}), event_list = {event_in};
    else, event_list = event_in; end

    % CPU timing accumulator
    cpu_total = 0;  n_secondary = 0;

    for k = 1:steps
        t = (k-1)*dt;

        % Events
        while ev <= numel(event_list) && t >= event_list{ev}{1}
            kind = event_list{ev}{2}; pl = event_list{ev}{3};
            if strcmp(kind,'load'), Pload(pl(1)) = Pload(pl(1)) + pl(2); end
            ev = ev + 1;
        end

        % RK4 on fast plant
        inj = @(uu,dwl) max(min(uu - p.KDR.*dwl, Pmax), Pmin) .* gen;  % rebind each step (stale-closure hygiene; load-only events here)
        f = @(de,dwl) fastderiv(de,dwl,u,gen,bus,lines,Pload,Pmax,p,inj);
        [a1,b1]=f(delta,dw);
        [a2,b2]=f(delta+dt/2*a1,dw+dt/2*b1);
        [a3,b3]=f(delta+dt/2*a2,dw+dt/2*b2);
        [a4,b4]=f(delta+dt*a3,  dw+dt*b3);
        delta = delta + dt/6*(a1+2*a2+2*a3+a4);
        dw    = dw    + dt/6*(b1+2*b2+2*b3+b4);

        % Slow secondary (every nsub steps)
        if mod(k-1, nsub) == 0
            tc = tic;
            Pi   = inj(u, dw);
            r    = zeros(1,N);
            r(gen) = Pi(gen) ./ max(Pmax(gen), 1e-9);
            du   = zeros(1,N);

            switch method

                case 'proposed'
                    % ------------------------------------------------
                    % PROPOSED: harmonic-mean Perron consensus on r_i
                    % du = ALPHA*(omega*-omega_i)*TS + KSHARE*(W*r - r)*Pmax
                    % (identical to paper1_reproduce.m lines 108-110)
                    % ------------------------------------------------
                    rmix = (W_harm * r')';
                    du(gen) = p.ALPHA*(-dw(gen))*p.TS ...  % local dw_i per Eq. (2c)
                            + p.KSHARE*(rmix(gen)-r(gen)).*Pmax(gen);

                case 'simpsonporco'
                    % ------------------------------------------------
                    % SIMPSON-PORCO et al. 2015 [simpson2] - frequency DAPI,
                    % their eq. (6):
                    %   k_I * dOmega_i/dt = -(omega_i - omega*)
                    %                       - sum_j a_ij*(Omega_i - Omega_j)
                    %   droop: omega_i = omega* - m_i P_i + Omega_i,  m_i = 1/KDR_i,
                    % which maps to this framework's setpoint as
                    %   u_i = KDR_i * Omega_i.
                    % At equilibrium: omega = omega* (exact restoration) and
                    % Omega in consensus, so P_i is proportional to 1/m_i = KDR_i
                    % (their Theorem): sharing is DROOP-proportional, and
                    % rating-proportional ONLY if droops are co-designed with
                    % m_i ~ 1/Pmax_i. The testbed droops are not, so a per-unit
                    % loading spread remains and the small-rating unit is driven
                    % to its limit -- exactly the design coupling that per-unit
                    % consensus (Theorem 3) removes. Metropolis weights.
                    % ------------------------------------------------
                    kI      = 1.0;                       % their integral gain
                    L_metro = diag(sum(A_metro,2)) - A_metro;
                    dOmega  = zeros(1,N);
                    dOmega(gen) = (1/kI) * ...
                        (-dw(gen) - (L_metro(gen,gen)*Omega(gen)')');
                    Omega(gen) = Omega(gen) + p.TS * dOmega(gen);
                    u_dapi  = p.KDR .* Omega;
                    du(gen) = u_dapi(gen) - u(gen);

                case 'ding'
                    % ------------------------------------------------
                    % DING et al. 2019 [ding-event]
                    % Their eq. (3)-(4) (time-triggered version for fair comparison):
                    %   u_dot_i = c1*(omega* - omega_i) + c2*sum_j a_ij*(r_j - r_i)
                    %
                    % Same consensus variable as proposed: r_i = P_i/P_i^max
                    % KEY DIFFERENCE: unweighted graph (a_ij = 1 on comm edges)
                    % => lower lambda2 => slower convergence
                    %
                    % Discretised at Ts = 0.02s (same as proposed)
                    % Gains c1, c2 matched to proposed ALPHA, KSHARE*eps for fair comparison
                    % ------------------------------------------------
                    rmix_uw = (W_unw * r')';
                    du(gen) = p.ALPHA*(-dw(gen))*p.TS ...  % local dw_i, per their eq. (3)-(4)
                            + p.KSHARE*(rmix_uw(gen)-r(gen)).*Pmax(gen);
                    % Note: same formula as proposed, only W_unw differs

                case 'shi'
                    % ------------------------------------------------
                    % SHI et al. 2020 [shi-piconsensus]
                    % Their frequency controller (eq. 7-8):
                    %   u_dot_i = kP*(omega* - omega_i)
                    %           + kI * xi_i
                    %           + c * sum_j a_ij*(r_j - r_i)
                    %   xi_dot_i = omega* - omega_i   (integral of freq error)
                    %
                    % Same consensus variable as proposed: r_i = P_i/P_i^max
                    % KEY DIFFERENCE: extra integral term + Metropolis weights
                    % => integral improves freq accuracy but adds dynamics
                    % => Metropolis weights => lower lambda2 than harmonic
                    % ------------------------------------------------
                    xi = xi + p.TS*(-dw);   % xi_dot = omega* - omega_i = -dw_i
                    rmix_shi = (W_metro * r')';
                    kP = p.ALPHA; kI = 0.05; c = p.KSHARE;
                    du(gen) = kP*(-dw(gen))*p.TS ...  % local dw_i, per their eq. (7)-(8)
                            + kI*xi(gen)*p.TS ...
                            + c*(rmix_shi(gen)-r(gen)).*Pmax(gen);
            end

            % Anti-windup + saturation (same for all)
            hi = (Pi >= Pmax-1e-6) & (du > 0);
            lo = (Pi <= Pmin+1e-6) & (du < 0);
            du(hi) = 0; du(lo) = 0;
            u = max(min(u + du, Pmax), Pmin);

            cpu_total = cpu_total + toc(tc)*1e6;  % microseconds
            n_secondary = n_secondary + 1;
        end

        tvec(k) = t;
        Fmat(k,:) = p.f0 + dw/(2*pi);
        Pi_s = inj(u,dw);
        rr = nan(1,N);
        rr(gen) = Pi_s(gen) ./ max(Pmax(gen),1e-9);
        Rmat(k,:) = rr;
    end

    R.t       = tvec;
    R.F       = Fmat;
    R.R       = Rmat;
    R.gen     = gen;
    R.lambda2 = l2;
    R.cpu_us  = cpu_total / max(n_secondary,1);
end


%% ========================================================================
function make_figure(RR, labels, lam2v, settle, niter, p, t_step) %#ok<INUSD>
    % Two-panel version for publication: (a) COI frequency, (b) sharing error.
    % Settling bars removed (Sec. 10.9 carries the update counts); DAPI's
    % non-settling is annotated textually instead of with sentinel bars.
    cols = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556];
    lsty = {'-','-','--',':'};
    lw   = [2.2 1.5 1.5 1.5];
    gen_idx = find(p.gen0);

    figure('visible','off','position',[100 100 700 520]);

    % Panel (a): COI frequency
    subplot(2,1,1); hold on; box on; grid on;
    for m2 = 1:4
        fCOI = mean(RR{m2}.F(:,gen_idx), 2);
        plot(RR{m2}.t, fCOI, lsty{m2}, 'Color',cols(m2,:), 'LineWidth',lw(m2), ...
             'DisplayName', labels{m2});
    end
    xline(t_step,'--','Color',[.5 .5 .5],'HandleVisibility','off');
    yline(50,':','Color',[.4 .4 .4],'HandleVisibility','off');
    ylabel('COI frequency (Hz)'); ylim([49.94 50.03]); xlim([0 40]);
    title('(a) Frequency restoration: all four methods reach 50 Hz');
    legend('Location','southeast','FontSize',7);

    % Panel (b): per-unit sharing error
    subplot(2,1,2); hold on; box on; grid on;
    for m2 = 1:4
        se = max(RR{m2}.R,[],2,'omitnan') - min(RR{m2}.R,[],2,'omitnan');
        plot(RR{m2}.t, se, lsty{m2}, 'Color',cols(m2,:), 'LineWidth',lw(m2), ...
             'HandleVisibility','off');
    end
    xline(t_step,'--','Color',[.5 .5 .5],'HandleVisibility','off');
    yline(0.01,':','Color',[.5 .5 .5],'HandleVisibility','off');
    text(36, 0.03, '1% threshold','FontSize',7,'Color',[.4 .4 .4]);
    text(20, 0.24, 'DAPI: droop-proportional equilibrium', ...
         'FontSize',8,'Color',cols(2,:));
    text(20, 0.205, '(0.28 spread; smallest-rated unit at its rating)', ...
         'FontSize',7,'Color',cols(2,:));
    text(4.3, 0.06, 'r_i-consensus methods: <1% in 9 updates', ...
         'FontSize',7,'Color',[.3 .3 .3]);
    ylabel('sharing error  max r_i - min r_i'); xlabel('time (s)');
    ylim([0 0.45]); xlim([0 40]);
    title('(b) Proportional sharing: consensus variable decides the equilibrium');

    saveas(gcf,'fig_comparison_correct.png');
end


%% ========================================================================
function write_table(labels, cites, lam2v, nadir, rocof, settle, niter, ...
                     sharerr, cpu_us)
    fid = fopen('tab_comparison_correct.tex','w');
    fprintf(fid,'\\begin{table*}[!tbp]\n');
    fprintf(fid,['\\caption{Quantitative comparison of four distributed secondary-control\n'...
        'methods that solve the same problem (frequency restoration to nominal\n'...
        'and proportional active-power sharing in islanded AC microgrids) under\n'...
        'Scenario~1 conditions ($+30$~MW load step at $t=2$~s, six-unit test system).\n'...
        'All methods share the same plant, primary droop gains, and communication\n'...
        'graph $\\mathcal{G}_C$; only the secondary control law and edge-weighting differ.\n'...
        '$\\lambda_2(L_C)$: algebraic connectivity (governs sharing convergence rate).\n'...
        'Settling: first instant the per-unit sharing error drops below $1\\%%$.\n'...
        'Updates: number of secondary control periods ($T_s=20$~ms) from disturbance\n'...
        'to settling. CPU: average per-update computation time.}\n']);
    fprintf(fid,'\\label{tab:comparison}\n\\centering\n');
    fprintf(fid,'\\setlength{\\tabcolsep}{3.5pt}\\footnotesize\n');
    fprintf(fid,'\\begin{tabular}{lcccccccc}\n\\toprule\n');
    fprintf(fid,['Method & Ref. & $\\lambda_2$ & Nadir & RoCoF & Settling & Updates & '...
                 'Final sharing & CPU\\\\\n']);
    fprintf(fid,'       &      & $(L_C)$ & (Hz) & (Hz/s) & (s) & to settle & error & ($\\mu$s)\\\\\n');
    fprintf(fid,'\\midrule\n');
    for m = 1:4
        if m==1, bo='\\textbf{'; bc='}'; else, bo=''; bc=''; end
        if isnan(settle(m)), st_s='$>40$'; else
            st_s=sprintf('%s%.1f%s',bo,settle(m),bc); end
        if isnan(niter(m)), ni_s='$>2000$'; else
            ni_s=sprintf('%s%d%s',bo,niter(m),bc); end
        fprintf(fid,'%s%s%s & %s & %s%.4f%s & %s%.4f%s & %s%.4f%s & %s & %s & %s%.2e%s & %s%.1f%s\\\\\n',...
            bo,labels{m},bc, cites{m},...
            bo,lam2v(m),bc, bo,nadir(m),bc, bo,rocof(m),bc,...
            st_s, ni_s,...
            bo,sharerr(m),bc, bo,cpu_us(m),bc);
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table*}\n');
    fclose(fid);
end


%% ========================================================================
%  Shared infrastructure — identical to paper1_reproduce.m
%% ========================================================================
function p = params()
    p.f0=50; p.w0=2*pi*p.f0; p.N=6;
    p.H=[5 4 6 4.5 5.5 4]; p.M=2*p.H/p.w0;
    p.Dvir=1.2*ones(1,p.N);
    p.Pmax=[.30 .25 .20 .20 .25 .20]; p.Pmin=zeros(1,p.N);
    p.KDR=[.20 .18 .20 .18 .20 .18];
    p.Pload0=[.08 .16 .12 .12 .14 .10];
    p.ALPHA=0.30; p.KSHARE=0.8; p.TS=0.02; p.RX=0.3;
    X=1/8; R=p.RX*X; den=R^2+X^2; p.G_L=R/den; p.B_L=X/den;
    p.LINES=[1 2;2 3;3 4;4 5;5 6;6 1;1 3;2 5];
    p.COMM =[1 2;2 3;3 4;4 5;5 6;6 1;1 3;1 4];
    p.gen0=[true true true true true false];
end
function w = pmean(a,b,pp)
    if pp==0, w=sqrt(a*b); else, w=((a^pp+b^pp)/2)^(1/pp); end
end
function A = Amat(edges,act,scheme,N)
    d=zeros(1,N);
    for e=1:size(edges,1)
        i=edges(e,1); j=edges(e,2);
        if act(i)&&act(j), d(i)=d(i)+1; d(j)=d(j)+1; end
    end
    A=zeros(N);
    for e=1:size(edges,1)
        i=edges(e,1); j=edges(e,2);
        if ~(act(i)&&act(j)), continue; end
        di=d(i); dj=d(j);
        switch scheme
            case 'harmonic',  w=1/pmean(di,dj,-1);
            case 'metropolis',w=1/(max(di,dj)+1);
            case 'localdeg',  w=1/max(di,dj);
        end
        A(i,j)=w; A(j,i)=w;
    end
end
function A = Amat_ones(edges,act,N)
    A=zeros(N);
    for e=1:size(edges,1)
        i=edges(e,1); j=edges(e,2);
        if act(i)&&act(j), A(i,j)=1; A(j,i)=1; end
    end
end
function W = Wperron(A,act)
    N=size(A,1); L=diag(sum(A,2))-A; idx=find(act);
    if numel(idx)<2, W=eye(N); return; end
    ep=0.6/max(max(sum(A(idx,:),2)),1e-9); W=eye(N)-ep*L;
end
function l2 = lam2_val(A,act)
    idx=find(act);
    if numel(idx)<2, l2=0; return; end
    As=A(idx,idx); L=diag(sum(As,2))-As; ev=sort(eig(L)); l2=ev(2);
end
function [dd,ddw]=fastderiv(delta,dw,u,gen,bus,lines,Pload,Pmax,p,inj)
    N=p.N; Pi=inj(u,dw); pe=zeros(1,N);
    for e=1:size(lines,1)
        i=lines(e,1); j=lines(e,2);
        if ~(bus(i)&&bus(j)), continue; end
        th=delta(i)-delta(j);
        pe(i)=pe(i)+p.G_L*(1-cos(th))+p.B_L*sin(th);
        pe(j)=pe(j)+p.G_L*(1-cos(-th))+p.B_L*sin(-th);
    end
    wbar=mean(dw(bus)); damp=p.Dvir.*(dw-wbar).*gen;
    dd=dw.*bus; ddw=((Pi-Pload-pe-damp)./p.M).*bus;
end