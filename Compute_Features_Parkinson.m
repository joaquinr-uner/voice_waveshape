root = '/media/Datos/joaquinruiz/PC-GITA/Matfiles/Results';

groups = {'Patologicas';'Control'};
vowels = {'A';'E'};
nv = 20;
r_opt = 15;outn = 1;
nmn = 1;

Nm = 44100;
for g=1:size(groups,1)
    group = groups{g};
    for v=1:size(vowels,1)

        vow = vowels{v};

        load(fullfile(root,['Results_' group '_' vow '_tvwse.mat']));

        J = length(Sis.Error_tvwse);

        dalpk = zeros(J,r_opt-1);
        s = zeros(J,r_opt);
        TV_B1 = zeros(1,J);
        TV_F0 = zeros(1,J);
        Error = Sis.Error_tvwse;
        for j=1:J
            v_ie = Sis.Coefs_tvwse(:,j);
            B1 = Sis.Amp(:,j);
            phi = Sis.Phi(:,j);
            N = length(B1);
            Np = floor(0.2*N);
            Next = N + 2*Np;
            ft = medfilt1(diff(phi),100);
            ft(end+1) = ft(end);
            %B1 = v1*A;
            [ti, alp, gamh, eh] = parse_coefs(Next,v_ie,r_opt,nv,0,outn);
            [ti,alp] = remove_outn(ti,alp,outn);
            [a, b, Q, q] = compute_hafs(ti,alp,gamh,'pchip',nmn,B1);
            %figure(1)
            %subplot(211)
            %plot(Sis.Signals(l,:))
            %subplot(212)
            %plot(Q')
            %title(['Signal ' num2str(Sis.Indx(l)) '. Type ' num2str(St.Types(l))])
            %pause(0.5)
            for m=1:r_opt-1
                tij = ti{m};
                difftk(m) = sum(abs(diff(tij(2:end-1))))/N;
                alpm = interp1(tij,q{m},1:length(B1),'pchip');
                dalpk(j,m) = sum(abs(diff(alpm)));
            end
            dB1 = diff(B1);
            dB1 = [dB1; dB1(end)];
            dB1 = medfilt1(dB1,100);
            TV_B1(j) = sum(abs(dB1))/length(dB1);
            TV_F0(j) = sum(abs(diff(ft)))/length(ft);
            s(j,:) = sqrt(a.^2+b.^2);
        end
        Feats = struct('Error_tvWSE',Error,'TV_F0',TV_F0,'TV_B1',TV_B1,'TVAmpl_k',dalpk,'s',s,'Ropt',...
            r_opt);
        save(fullfile(root,['Features_' group '_' vow '.mat']),'Feats')

        TabFeats= table(Error',TV_B1',TV_F0',dalpk,s);
        TabFeats.Properties.VariableNames = {'Error','TV_B1','TV_F0','TV_HAF','S'};
        writetable(TabFeats,['Features_' group '_' vow '.csv'])
    end
end