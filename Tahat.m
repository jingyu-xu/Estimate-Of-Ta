fid = fopen('breastinfo_simple.txt','r');
breastID = str2double(fgetl(fid));
s1 = str2double(fgetl(fid));
s2 = str2double(fgetl(fid));
s3 = str2double(fgetl(fid));
class = str2double(fgetl(fid));
fclose(fid);

load mtype.mat;
load pval.mat;

muscle_wall = 153;
skin_start = 138;

% Convert vector into cube
mtype_cube = zeros(s1,s2,s3); % each voxel is .5mmx.5mmx.5mm
pval_cube = zeros(s1,s2,s3);
cur_pos = 1;
for k=1:s3
    for j=1:s2
        for i= 1:s1
            mtype_cube(i,j,k) = mtype(cur_pos);
            pval_cube(i,j,k) = pval(cur_pos);
            cur_pos = cur_pos + 1;
        end 
    end
end

% subsample cubes in order to solve sparse matrix
s1_ss = floor(s1/2); % voxels are now 1mmx1mmx1mm
s2_ss = floor(s2/2);
s3_ss = floor(s3/2);
xi = 1; yi = 1; zi = 1;
mtype_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
pval_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
for z=1:2:s3-1
    for y = 1:2:s2
        for x = 1:2:s1
            mid = mtype_cube(x,y,z);
            pid = pval_cube(x,y,z);
            mtype_cube_subsamp(xi,yi,zi) = mid;
            pval_cube_subsamp(xi,yi,zi) = pid;
            xi = xi+1;
        end
        xi = 1;
        yi = yi + 1;
    end
    yi = 1;
    zi = zi + 1;
end
% some voxels not converted to muscle during subsampling
% so do that now
% still need to figure out how to get pval converted for fdtd
for x=1:s1_ss
    for y=1:s2_ss
        for z=1:s3_ss
            if x > 153
               mtype_cube_subsamp(x,y,z) = -4;
                pval_cube_subsamp(x,y,z) = 1;
            end
        end
    end
end

% save mtype_cube_subsamp; save pval_cube_subsamp.mat;
% load mtype_cube_subsamp.mat;
[s1_ss, s2_ss, s3_ss] = size(mtype_cube_subsamp);
model = mtype_cube_subsamp; air_id = -1;
% sliceomatic(model)
% model = model(18:end,:,:);
% [s1_ss, s2_ss, s3_ss] = size(model);

tumor_on = 0; tumor_depth = 8;  tumor_radius = 10; Tambient = 27; Tart = 37; 
[T_3d_nom,tissue_3d_nom] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start);
figure; contourf(T_3d_nom(:,:,floor(s3_ss/2)));
colorbar;
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
title('Normal Temperature Profile (\circC)','FontSize',14);

% Generate temperature anomalies with radius = 10
tumor_on = 1; tum_y_cen = 90; tum_z_cen = floor(s3_ss/2); 
tumor_radius = 10; tumor_depth = 10; tum_x_cen = 20 + tumor_depth;% tumor dept is 1cm
[T_3d_abn1,tissue_3d_abn1] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);
tumor_radius = 10; tumor_depth = 20; tum_x_cen = 20 + tumor_depth;% tumor depth is 2cm
[T_3d_abn2,tissue_3d_abn2] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);
tumor_radius = 10; tumor_depth = 30; tum_x_cen = 20 + tumor_depth;% tumor depth is 3cm
[T_3d_abn3,tissue_3d_abn3] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);
tumor_radius = 10; tumor_depth = 40; tum_x_cen = 20 + tumor_depth;% tumor depth is 4cm
[T_3d_abn4,tissue_3d_abn4] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);

% Plot temperature change due to tumor
T_diff1 = T_3d_abn1 - T_3d_nom;
T_diff2 = T_3d_abn2 - T_3d_nom;
T_diff3 = T_3d_abn3 - T_3d_nom;
T_diff4 = T_3d_abn4 - T_3d_nom;
plot(T_diff1(20:100,90,80),'k','LineWidth',1);
load scat_fib_nom_22ghz_ceramic_center1.mat;
load scat_fib_nom_22ghz_ceramic_right1.mat;
load scat_fib_nom_22ghz_ceramic_left1.mat;
load scat_fib_nom_22ghz_ceramic_top1.mat;
load scat_fib_nom_22ghz_ceramic_bottom1.mat;

model = model(18:end,:,:);
[s1_ss, s2_ss, s3_ss] = size(model);

% Read WFs into a cube
wf_cube_2_2ghz_ceramic_center = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_center1);
wf_cube_2_2ghz_ceramic_right = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_right1);
wf_cube_2_2ghz_ceramic_left = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_left1);
wf_cube_2_2ghz_ceramic_top = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_top1);
wf_cube_2_2ghz_ceramic_bottom = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_bottom1);

% eliminate extra data
wf_cube_2_2ghz_ceramic_center = wf_cube_2_2ghz_ceramic_center(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_right = wf_cube_2_2ghz_ceramic_right(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_left = wf_cube_2_2ghz_ceramic_left(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_top = wf_cube_2_2ghz_ceramic_top(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_bottom = wf_cube_2_2ghz_ceramic_bottom(:,1:s2_ss,:);

% plot 3d
% sliceomatic(wf_cube_2_2ghz_ceramic_center);
% find norms of wfs and plot cube results
% convert cube to 1d
wf_vec_2_2ghz_ceramic_center = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_center,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_right = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_right,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_left = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_left,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_top = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_top,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_bottom = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_bottom,s1_ss,s2_ss,s3_ss);

% sum and find the norm of each wf
wf_vec_2_2ghz_norm_ceramic_center = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_center);
wf_vec_2_2ghz_norm_ceramic_right = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_right);
wf_vec_2_2ghz_norm_ceramic_left = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_left);
wf_vec_2_2ghz_norm_ceramic_top = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_top);
wf_vec_2_2ghz_norm_ceramic_bottom = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_bottom);

wf_cube_2_2ghz_norm_ceramic_center = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_center,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_right = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_right,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_left = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_left,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_top = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_top,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_bottom = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_bottom,s1_ss,s2_ss,s3_ss);
% construcate weighting function matrix
WF = [wf_vec_2_2ghz_norm_ceramic_center wf_vec_2_2ghz_norm_ceramic_right wf_vec_2_2ghz_norm_ceramic_left ...
    wf_vec_2_2ghz_norm_ceramic_top wf_vec_2_2ghz_norm_ceramic_bottom];

Tvec_nominal = convert_3d_to_1d(T_3d_nom,s1_ss,s2_ss,s3_ss);

TB_nominal = (WF'*Tvec_nominal);


% Calculate Brightness Temperature with Tumors
Tvec_abn1 = convert_3d_to_1d(T_3d_abn1,s1_ss,s2_ss,s3_ss);
Tvec_abn2 = convert_3d_to_1d(T_3d_abn2,s1_ss,s2_ss,s3_ss);
Tvec_abn3 = convert_3d_to_1d(T_3d_abn3,s1_ss,s2_ss,s3_ss);
Tvec_abn4 = convert_3d_to_1d(T_3d_abn4,s1_ss,s2_ss,s3_ss);
TB_abn1 = (WF'*Tvec_abn1);
TB_abn2 = (WF'*Tvec_abn2);
TB_abn3 = (WF'*Tvec_abn3);
TB_abn4 = (WF'*Tvec_abn4);

% Calculate Ta, Ta = actual temperature - nominal temperatur with tumor
Ta1 = Tvec_abn1 - Tvec_nominal;
Ta2 = Tvec_abn2 - Tvec_nominal;
Ta3 = Tvec_abn3 - Tvec_nominal;
Ta4 = Tvec_abn4 - Tvec_nominal;

% Calculate Ta0 without tumor
Ta0 = 0;

% Generate White Gaussian Noise TBn
std = 0.1;
elps = 1e-8;
TBn = random('norm',0,0.1,5,1);
figure
plot(TBn);

% Calculate eignvector and eignvalue of w'w and ww'
[theta, D] = eig(WF'*WF);
lambda = eig(WF'*WF);
lambda_R = eig(WF'*WF) + elps;
P = WF * theta;

% Calculate TB_delta TB delta = Brightness Temperature -  TB0 with tumor
% and without tumor TB_delta0
TB_delta0 = TBn;
TB_delta = TB_abn1 - TB_nominal;
TB_delta1 = TB_abn1 - TB_nominal + TBn;
TB_delta2 = TB_abn2 - TB_nominal + TBn;
TB_delta3 = TB_abn3 - TB_nominal + TBn;
TB_delta4 = TB_abn4 - TB_nominal + TBn;

% dot_abn without noise
dot11 = dot(TB_delta,theta(:,1));
dot12 = dot(TB_delta,theta(:,2));
dot13 = dot(TB_delta,theta(:,3));
dot14 = dot(TB_delta,theta(:,4));
dot15 = dot(TB_delta1,theta(:,5));
dot_vec = [dot11; dot12; dot13; dot14; dot15];
% dot_abn
dot1 = dot(TB_delta1,theta(:,1));
dot2 = dot(TB_delta1,theta(:,2));
dot3 = dot(TB_delta1,theta(:,3));
dot4 = dot(TB_delta1,theta(:,4));
dot5 = dot(TB_delta1,theta(:,5));
dot_vec_1 = [dot1; dot2; dot3; dot4; dot5];

% dot_nom
dot1_nom = dot(TB_delta0,theta(:,1));
dot2_nom = dot(TB_delta0,theta(:,1));
dot3_nom = dot(TB_delta0,theta(:,1));
dot4_nom = dot(TB_delta0,theta(:,1));
dot5_nom = dot(TB_delta0,theta(:,1));
dot_vec_nom = [dot1_nom; dot2_nom; dot3_nom; dot4_nom; dot5_nom];

lambda_inverse = 1./lambda;
lambda_R_inverse = 1./lambda_R;
Tahat_N = P * ((lambda_inverse).*dot_vec); % without tumor without regularize without noise
Tahat_R_N = P * ((lambda_R_inverse).*dot_vec); % without tumor with regularized without noise
Tahat_nom = P * ((lambda_inverse).*dot_vec_nom); % without tumor without regularize
Tahat_R_nom = P * ((lambda_R_inverse).*dot_vec_nom); % without tumor with regularized
Tahat_abn1 = P * ((lambda_inverse).*dot_vec_1); % with tumor without regularize
Tahat_R_abn1 = P * ((lambda_R_inverse).*dot_vec_1); % with tumor with regularize
Tahat_N_3d = convert_1d_to_3d(Tahat_N,s1_ss,s2_ss,s3_ss);
Tahat_R_N_3d = convert_1d_to_3d(Tahat_R_N,s1_ss,s2_ss,s3_ss);
Tahat_3d_nom = convert_1d_to_3d(Tahat_nom,s1_ss,s2_ss,s3_ss);
Tahat_R_3d_nom = convert_1d_to_3d(Tahat_R_nom,s1_ss,s2_ss,s3_ss);
Tahat_3d_abn1 = convert_1d_to_3d(Tahat_abn1,s1_ss,s2_ss,s3_ss);
Tahat_R_3d_abn1 = convert_1d_to_3d(Tahat_R_abn1,s1_ss,s2_ss,s3_ss);

% Plot unregularized Ta (with tumor
% without noise)
figure % Plot Ta without R with tumor without noise
plot(Tahat_N_3d(20:100,90,floor(s3_ss/2)),'k','LineWidth',1.5);
xlabel('Depth (mm)');ylabel('Tahat');
figure % Plot Ta with R with tumor without noise
plot(Tahat_R_N_3d(20:100,90,floor(s3_ss/2)),'r','LineWidth',1.5);
xlabel('Depth (mm)');ylabel('Tahat');


figure; % Plot Ta_R without Noise with tumor
plot(Tahat_R_N_3d(20:100,90,floor(s3_ss/2)),'k','LineWidth',1.5);
xlabel('Depth (mm)');ylabel('Tahat');
figure; % Plot Ta_R with Noise with tumor
plot(Tahat_R_3d_abn1(20:100,90,floor(s3_ss/2)),'k','LineWidth',1.5);
xlabel('Depth (mm)');ylabel('Tahat');


% Plot the difference of Ta with and without tumor(without regularized with Noise)
figure; % Plot Ta without R with Noise without Tumor
plot(Tahat_3d_nom(20:100,90,floor(s3_ss/2)),'k','LineWidth',1.5); %Tahat without tumor without regularize
xlabel('Depth (mm)');ylabel('Tahat');
figure; % Plot Ta without R with Noise with Tumor
plot(Tahat_3d_abn1(20:100,90,floor(s3_ss/2)),'r','LineWidth',1.5);  %Tahat with Tumor without regularize
xlabel('Depth (mm)');ylabel('Tahat');



figure; % Plot Ta_R with Noise without tumor
plot(Tahat_R_3d_nom(20:100,90,floor(s3_ss/2)),'k','LineWidth',1.5);
xlabel('Depth (mm)');ylabel('Tahat');
figure; % Plot Ta_R with Noise with tumor
plot(Tahat_R_3d_abn1(20:100,90,floor(s3_ss/2)),'r','LineWidth',1.5);
xlabel('Depth (mm)');ylabel('Tahat');


% Export useful variables to Excel
% xlswrite('Tvec_nominal.xlsx',Tvec_nominal);
% xlswrite('TB_nominal.xlsx',TB_nominal);
% xlswrite('TB_abn1.xlsx',TB_abn1);
% xlswrite('Ta1.xlsx',Ta1);
% xlswrite('TBn.xlsx',TBn);
% xlswrite('theta.xlsx',theta);
% xlswrite('lambda.xlsx',lambda);
% xlswrite('lambda_R.xlsx',lambda_R);
xlswrite('P.xlsx',P(1:5));
xlswrite('TB_delta1.xlsx',TB_delta1);
xlswrite('TB_delta.xlsx',TB_delta);
% xlswrite('Tahat_R.xlsx',Tahat_R);
% xlswrite('Tahat_R_nom.xlsx',Tahat_R_nom);
% xlswrite('Tahat_R_abn1.xlsx',Tahat_R_abn1);



















