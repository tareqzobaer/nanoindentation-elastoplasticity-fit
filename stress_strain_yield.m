warning('off', 'all')
%%  indenter data

% E: elastic modulus of the material
% v: Poisson's ratio of the material
% E_i: elastic modulus of the indenter
% v_i: Poisson's ratio of the indenter
% R: sphere radius

E_i = 345e9; % in Pa
v_i = 0.29;
v = 0.18;
R = 158/2*1e3; %in nm

yield_point_offset = 0.0/100;
elastic_modulus_strain = 0.5/100;
%% read experiment data
test_file_names = {'sp_sample1_1.xlsx','sp_sample1_2.xlsx','sp_sample2.xlsx'...
    'sp_sample3.xlsx','sp_sample5_1.xlsx','sp_sample5_2.xlsx','sp_sample6_1.xlsx','sp_sample7.xlsx'};


if ~isfolder('ind_se_data')
    mkdir('ind_se_data');
end

if isfile('ind_se_data\var_limits.mat')
    load('ind_se_data\var_limits.mat','var_limits')
end

% 
lim_f = zeros(length(test_file_names),1);
lim_h = zeros(length(test_file_names),1);
lim_s = zeros(length(test_file_names),1);
lim_e = zeros(length(test_file_names),1);
for j=1:length(test_file_names)
    test_file_name = strcat('data\',test_file_names{j});

    [~, sheet_names] = xlsfinfo(test_file_name);
    results_data = readmatrix(test_file_name, 'Sheet', sheet_names{1}, 'UseExcel', 0);

    sheet_names = sheet_names(1:length(sheet_names));
    test_num = length(sheet_names);


    tests = 1:test_num;
    results_data = results_data(tests+1,:);

    mat_prop = zeros(test_num,4);


    for i=tests
        sheet = sheet_names{i};
        sheet_num = str2double(sheet(end-2:end));
        unload_E = results_data(sheet_num,2); % 2: unload_modulus

        loading_data = readtable(test_file_name, 'Sheet', sheet, 'UseExcel', 0);
        loading_data = loading_data(:,1);
        load_smt = ~cellfun(@isempty,loading_data.Segment);

        segments = find(load_smt);

        load_start= segments(1);
        load_end = segments(2)-1;

        columns = [2 3 13 16 17 6 25];
        % columns = [h_t P S A h_c E];

        data = readmatrix(test_file_name, 'Sheet', sheet, 'UseExcel', 0);

        % data = data(load_start:load_end,columns);
        data = data(:,columns);
        nan_data = any(isnan(data),2);
        data(nan_data,:) = []; %the row with end is being deleted

        imaginary_data = data(:,4)<0;
        data(imaginary_data,:) = []; %the row with end is being deleted
       
        % S: the contact stiffness using CSM (measured)
        % h_t: total penetration depth (directly measured)
        % P: load (directly measured)
        % h_c: contact depth% 
        % A: projected contact area

        h_t = data(:,1).*1e-9; % nm --> m
        P = data(:,2).*1e-3; % mN -->m
        S = data(:,3); % N/m
        E_IT = data(:,6).*1e9; % GPa --> Pa
        E_r = data(:,7).*1e9; % GPa --> Pa
        A = data(:,4).*1e-18; % nm2 -->m2


        a = (sqrt(1/pi*A));
        h_c = data(:,5).*1e-9; % nm --> m

        psi  = 0.9.*(h_t./h_c - 1);
        s_indt = psi.*P./(pi.*a.^2);

        t1 = s_indt./(1-v^2);
        t2 = 4./( 2.7.*P.*( h_t./h_c - 1 ) );
        t3 = ( 2.*(h_t-h_c) ./ ( 2.*h_c./(a.^2 + h_c.^2) ).^(1/3) ).^(3/2);
        t4 = (1-v_i^2)./E_i;

        e_indt = (t1.*(t2.*t3 - t4));
        
        Emat = (1-v^2)./( ( (4./3./P).*t3 ) - t4 );
        
        S_e_pre = [s_indt e_indt];

        img_idx1=find(imag(S_e_pre(:,1))~=0);
        img_idx2=find(imag(S_e_pre(:,2))~=0);
        img_idx = union(img_idx2,img_idx2);
        S_e_pre(img_idx,:) = [];       

        nan_data = any(isnan(S_e_pre),2);
        S_e_pre(nan_data,:) = [];
        
        inf_data = any(S_e_pre>1e50,2);
        S_e_pre(inf_data,:) = [];


        [~,se_eps] = sort(S_e_pre(:,2));
        
        S_e = zeros(size(S_e_pre));
        S_e(:,1) = S_e_pre(se_eps,1);
        S_e(:,2) = S_e_pre(se_eps,2);
        
        % manual data cleanup
        if j==7 && i==1
            S_e([end-1 end],:) = [];
        elseif j==6 && i==2
            S_e([51 54 55 62 64 66 67 68 69 70 71 72 73 78 81 82 83 85 88 90 91 92 94 95 96 97 100 101 103 104 105 107 110 112 116 117 127 129 131 135 141 142 146 147 149 154 155 161],:) = [];
        end

        fitdata = polyfit(S_e(:,2),S_e(:,1)/1e6,12);
        stress_smooth = polyval(fitdata,S_e(:,2));

        
        stress = S_e(:,1)/1e6;
        strain = S_e(:,2);

        % Calculate first derivative (slope)
        d1 = gradient(stress_smooth, strain);

        % Calculate second derivative
        d2 = gradient(d1, strain);

        % Find where second derivative crosses zero (from positive to negative)
        zero_crossings = find(d2(1:end-1) .* d2(2:end) < 0 & d1(1:end-1) > 0);
        valid_zero_crossings = strain(zero_crossings)>0.015;
        zero_crossings(valid_zero_crossings)=[];

        % Find the maximum slope among these points
        [max_slope, max_index] = max(d1(zero_crossings));

        % The strain and stress at the point of maximum slope
        zero_crossing_strain = strain(zero_crossings(max_index));
        zero_crossing_stress = stress_smooth(zero_crossings(max_index));

        max_slope_strain = strain(zero_crossings(max_index))+yield_point_offset;
        max_slope_stress = polyval(fitdata,max_slope_strain);
        
        lin_strain = zero_crossing_strain;
        lin_part = S_e(:,2)<=lin_strain & S_e(:,2)>=lin_strain-elastic_modulus_strain;

        E = polyfit(S_e(lin_part,2),S_e(lin_part,1),1);
        m_l = E(1);
        b_l = E(2);
        
        mat_prop(i,1) = sheet_num;
        mat_prop(i,2) = m_l/1e9;
        mat_prop(i,3) = max_slope_stress;
        mat_prop(i,4) = max_slope_strain;
        
        fontsize = 15;
        fig = figure(1);
        title(strcat('SAMPLE : ',test_file_names{j}(10:end-4),sheet),'Interpreter','none','FontSize',fontsize)
        subplot(1,2,1)
        plot(data(:,1),data(:,2), 'LineWidth',2,'Color','k')
        xlabel('Indentation depth (nm)','FontSize',fontsize)
        ylabel('Load on sample (mN)','FontSize',fontsize)
        if isfile('ind_se_data\var_limits.mat')
            xlim([0 1.025*var_limits(j,2)])
            ylim([0 1.025*var_limits(j,1)])
        end
        box off
        title('Force-displacement curve','FontSize',fontsize)
        set(gca,'fontsize', fontsize,'fontname','times')
        
        
        subplot(1,2,2)
        scatter(S_e(:,2),S_e(:,1)/1e6,'Marker','.','LineWidth',3,'MarkerEdgeColor','k')
        hold on
        plot(S_e(:,2),polyval(fitdata,S_e(:,2)),'LineWidth',2,'Color','r')
        hold on
        plot(S_e(lin_part,2),polyval(E,S_e(lin_part,2))/1e6, 'LineWidth', 2,'Color','b')
        hold on
        scatter(max_slope_strain,max_slope_stress,'LineWidth',3,'MarkerEdgeColor','cyan','MarkerFaceColor','k' )
        % plot(offset_strain_data,offset_stress_se/1e6, 'LineWidth', 2)
        hold off

        xlabel('Indentation strain','FontSize',fontsize)
        ylabel('Indentation stress (MPa)','FontSize',fontsize)
        if isfile('ind_se_data\var_limits.mat')
            xlim([0 1.025*var_limits(j,4)])
            ylim([0 1.025*var_limits(j,3)])
        end
        legend('Data','Fitted curve', 'Linear part','Yield point','FontSize', fontsize+1,'Location','southeast');
        title('Stress-strain curve','FontSize',fontsize)
        set(gca,'fontsize', fontsize,'fontname','times')
        % text(0.005,230,strcat('Linear modulus, $E$ = ',num2str(round(mat_prop(j,1),2)),' GPa'),'FontSize',15,'Interpreter','latex','Color', 'b')
        % text(0.005,210,strcat('Initial yield strain, $\varepsilon_0$= ',num2str(round(mat_prop(j,3),4))),'FontSize',15,'Interpreter','latex','Color', 'b')
        % text(0.005,190,strcat('Initial yield stress, $\sigma_0$ = ',num2str(round(mat_prop(j,2),2)),' MPa'),'FontSize',15,'Interpreter','latex','Color', 'b')
    %     text(0.02, 0.98, strcat('Elastic modulus, $E$ = ',num2str(round(mat_prop(i,2),2)),' GPa'), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    % 'HorizontalAlignment', 'left', 'FontSize', 15, 'Interpreter','latex','Color', 'b');
    %     text(0.02, 0.94, strcat('Initial yield stress, $\sigma_y^0$ = ',num2str(round(mat_prop(i,3),2)),' MPa'), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    % 'HorizontalAlignment', 'left', 'FontSize', 15, 'Interpreter','latex','Color', 'b');
    %     text(0.02, 0.90, strcat('Initial yield strain, $\varepsilon_y^0$= ',num2str(round(mat_prop(i,4),4))), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    % 'HorizontalAlignment', 'left', 'FontSize', 15, 'Interpreter','latex','Color', 'b');

        % text(0.015,100,,'FontSize',15,)
        % text(0.015,80,strcat(sheet),'FontSize',15)
        % text(0.015,60,strcat('Yield calculation: ',num2str(offset_strain*100),'% offset method'),'FontSize',12,'Interpreter','none')
        % 
        % % text(0.001,190,strcat('$E_{unload}$ = ',num2str(unload_E),' GPa'),'FontSize',15,'Interpreter','latex')
        % % text(0.001,170,strcat('Initial yield stress, $\sigma_0$ = ',num2str(round(S_e(yield_data_unload,1)/1e6,2)),' MPa'),'FontSize',15,'Interpreter','latex')
        % % text(0.001,150,strcat('Initial yield strain, $\varepsilon_0$= ',num2str(S_e(yield_data_unload,2))),'FontSize',15,'Interpreter','latex')
        % % 
        % % % 
        % text(0.025,190,strcat('$E_{\sigma-\varepsilon}$ = ',num2str(m_l/1e9),' GPa'),'FontSize',15,'Interpreter','latex','Color', 'b')
        % text(0.025,170,strcat('Linear approximation till $\varepsilon$ = ',num2str(lin_strain)),'FontSize',15,'Interpreter','latex','Color', 'b')
        lim_f(j) = max(lim_f(j),max(P)*1e3);
        lim_h(j) = max(lim_h(j),max(h_t)*1e9);
        lim_s(j) = max(lim_s(j),max(S_e(:,1))/1e6);
        lim_e(j) = max(lim_e(j),max(S_e(:,2)));

        if ~isfolder('ind_se_data\figure')
            mkdir('ind_se_data\figure');
        end
        
        exportgraphics(fig,strcat('ind_se_data\figure\',test_file_names{j}(1:end-5),'_',sheet,'.png'),'Resolution',300)

       
        % pause

        se_table = table(strain,stress);
        writetable(se_table, strcat('ind_se_data\se_experiment',test_file_names{j}(3:end-5),'.xlsx'), 'WriteVariableNames', true,'Sheet',sheet);
    end
    
     writematrix( flipud(mat_prop),strcat('ind_se_data\sp_ind_elastic',test_file_names{j}(3:end-5),'.csv'))
    % hold off
    % saveas(fig,strcat(test_file_names{j}(10:end-5),'_',sheet,'.png'))
    % disp(strcat('E__',test_file_names{j}(1:end-4),' : ',num2str(mean(average_E)/1e9), char(177),num2str(std(average_E)/1e9),' GPa'));
    % % pause
end

var_limits = [lim_f lim_h lim_s lim_e];
save('ind_se_data\var_limits.mat','var_limits')




