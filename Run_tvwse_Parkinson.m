root = '/media/Datos/joaquinruiz/PC-GITA/Matfiles/Vowels';

groups = {'Patologicas';'Control'};
vowels = {'A';'E'};
mInterp = 'pchip';
nv = 20;
r_opt = 15;

sigma = 3e-6;
fmax = 0.05;
Crit = 'Wang';
CritParams = [2.1,4,6,12];
Nm = 44100;
for g=1:size(groups,1)
    group = groups{g};
    for v=1:size(vowels,1)

        vow = vowels{v};

        drt = dir(fullfile(root,group,vow));

        drt = drt(3:end);

        J = length(drt);

        Signals = zeros(Nm,J);
        S_tvWSE = zeros(Nm,J);
        Amp = zeros(Nm,J);
        Phi = zeros(Nm,J);
        Len = zeros(1,J);
        Coefs = zeros((r_opt-1)*(2*nv+4),J);
        E_tvwse = zeros(1,J);
        T_tvwse = zeros(1,J);
        for j=1:J
            name = drt(j).name;
            load(fullfile(root,group,vow,name))
            fprintf(['Running algorithm on signal ' name '...\n'])

            st = 6e3;
            N = min([length(data)-st,fs]);
            b = round(3/pi*sqrt(sigma/2)*N);
            s = data(st:st+N-1)';
            [F, sF] = STFT_Gauss(s,N,sigma,fmax);

            t = 0:1/fs:N/fs - 1/fs;
            f = 0:fs/N:fs*fmax - fs/N;
            U = istct_fast(F,f,0.3);

            W = U.*F;

            c = ridge_ext(W,0.1,0.1,10,10);

            [A_est, phi_est] = extract_harmonics(F,sF,c,b,b,1);

            fprintf(['Processing ' vow ' ' group '. Signal ' num2str(j) '...\n'])

            r_max = floor((fs/2)/f(max(c)));
            bwM = r_max*r_opt + 200*1/N;

            C = construct_dct(A_est,phi_est,r_opt);
            v = ((C'*C)\C')*s';

            s_norm = s./(v(1)*A_est);
            if bwM <fmax
                fmi = fmax;
            else
                if bwM > 0.5
                    fmi = 0.5;

                else
                    fmi = bwM;
                end
            end

            [Fn,sFn] = STFT_Gauss(s_norm,N,sigma,fmi);

            vnv = NNodes(nv,r_opt);

            Np = 0.1*N;

            Next = N + 2*Np;

            be = round(3/pi*sqrt(sigma/2)*Next);
            s_ext = extendSig(s,phi_est,3,Np,'fw-bw');

            s_ext = s_ext - mean(s_ext);

            [Fext, sFext] = STFT_Gauss(s_ext,Next,sigma,fmax);

            t_ext = 0:1/fs:Next/fs-1/fs;
            f_ext = 0:fs/Next:fs*fmax-fs/Next;
            Uext = istct_fast(Fext,f_ext,0.3);

            Wext = Uext.*Fext;

            cext = ridge_ext(Wext,0.1,0.1,10,10);

            cext = ridge_correct(cext,Fext,b,1);

            [A_ext,phi_ext] = extract_harmonics(Fext,sFext,cext,be,be,1);

            Ce = construct_dct(A_ext,phi_ext,r_opt);
            ve = ((Ce'*Ce)\Ce')*s_ext;

            se_norm = s_ext./(ve(1)*A_ext');

            Cen = construct_dct(ones(1,Next),phi_ext,r_opt);
            ven = ((Cen'*Cen)\Cen')*se_norm;

            vh = Init_tvWSE(ven,vnv,r_opt,1,N,Next);

            [lb,ub] = create_bounds(vnv,vh,r_opt,Next,1);

            tic;
            [sen_tvwse_n,v_ie,eflag_tvwse] = tvWSE(se_norm',ones(1,Next),phi_ext,r_opt,vnv,vh,mInterp,lb,ub,1,1,1);

            se_tvwse = ve(1)*A_ext.*sen_tvwse_n;
            if ~isrow(se_tvwse)
                se_tvwse = se_tvwse';
            end
            t_tvwse = toc;

            s_tvwse = se_tvwse(Np+1:N+Np);
            A = A_ext(Np+1:N+Np);
            phi = phi_ext(Np+1:N+Np);

            if N<Nm
                s = [s; zeros(Nm-N,1)];
                s_tvwse = [s_tvwse; zeros(Nm-N,1)];
                A = [A; zeros(Nm-N,1)];
                phi = [phi; zeros(Nm-N,1)];
            end

            e_tvwse = norm(s-s_tvwse)/norm(s);

            Signals(:,j) = s;
            S_tvWSE(:,j) = s_tvwse;
            Amp(:,j) = A;
            Phi(:,j) = phi;
            Coefs(:,j) = v_ie;
            Len(j) = N;
            E_tvwse(j) = e_tvwse;
            T_tvwse(j) = t_tvwse;
        end
        Sis = struct('Signals',Signals','Amp',Amp,...
                        'Phi',Phi,'S_tvwse',S_tvWSE,'Error_tvwse',...
                        E_tvwse,'Coefs_tvwse',Coefs,...
                        'Times_tvwse',T_tvwse,'Len',Len,'r_fix',r_opt,...
                        'I_fix',nv);
        save(fullfile(drt_r,['Results_' group '_' vow '_tvwse.mat']),'Sis')
    end
end
