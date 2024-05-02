addpath(genpath('/home/joaquinruiz/Dropbox/Github'))
addpath(genpath('C:\Users\Intel\Dropbox\Github'))
addpath(genpath('/home/jruiz'))
drt_r = 'E:\Voces\Results\VoiceTyping';
%drt_r = '/media/sentey/TOSHIBA EXT/Voces/Results/VoiceTyping';
%drt_r = '/home/joaquinruiz/Documents/Experimentos Saar Meei/Results_tvwse_Jun22';
%drt_r = 'Results';

mInterp = 'pchip';

DB = 'Meei';
dw = 0.04;
sigma = 1*(dw/(3))^2/64;
fmax = 0.05;

load('datos_saar.mat')

types = datos_saar.types;
fs = datos_saar.fs;
data = datos_saar.signals;

[Mm,N] = size(data);

b = round(3/pi*sqrt(sigma/2)*N);
t = 0:1/fs:N/fs - 1/fs;
f = 0:fs/N:fs*fmax - fs/N;

r_opt = 18;

nv = 15;

Signals = zeros(N,Mm);
S_tvWSE = zeros(N,Mm);
Amp = zeros(N,Mm);
Phi = zeros(N,Nm);
Coefs = zeros((r_opt-1)*(2*nv+4),Mm);
E_tvwse = zeros(1,Mm);
T_tvwse = zeros(1,Mm);
for i=1:Mm
    s = data(i,:);
    s = s - mean(s);
    type = types(i);

    [F, sF] = STFT_Gauss(s,N,sigma,fmax);

    U = istct_fast(F,f,0.3);

    W = U.*F;

    c = ridge_ext(W,0.1,0.1,10,10);

    [A_est, phi_est] = extract_harmonics(F,sF,c,b,b,1);

    fprintf(['Processing ' DB ' Signal ' num2str(i) '...\n'])

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

    e_tvwse = norm(s-s_tvwse)/norm(s);

    Signals(:,i) = s;
    Amp(:,i) = A;
    Phi(:,i) = phi;
    S_tvWSE(:,i) = s_tvwse;
    Coefs(:,i) = v_ie;
    E_tvwse(i) = e_tvwse;
    T_tvwse(i) = t_tvwse;

    fprintf(['Processing Completed for Voice Signal ' i '. Type: ' num2str(type) '. tvWSE Completed. Time: ' num2str(t_tvwse) '\n'...
        '. tvWSE Estimation Error: ' num2str(e_tvwse) '. lsqcurvefit exit flag: ' num2str(eflag_tvwse) '\n'])
end
Sis = struct('Signals',Signals','Types',Types,'Amp',Amp','Phi',Phi,...
    'S_tvwse',S_tvwse,'Error_tvwse',E_tvwse,'Coefs_tvwse',Coefs,...
    'Times_tvwse',T_tvwse,'r_fix',r_opt,'I_fix',nv);
save(fullfile(drt_r,['Results_tvwse_' DB '.mat']),'Sis')