% ʹ�ó�����ܣ�
% 1���޸���Ʊ���������NPar�������������Ʊ���������
% 2���޸���Ʊ����������ޣ�VarLow ��VarHign���������Ʊ����������ޣ�������ά������1����������һ�¡�
% 3���޸�FunName��������Ż�Ŀ�꺯��ֵ�ļ��㺯����
% 4���޸�����������MaxIterationsͨ�������ҳ����ʵ�MaxIterations����


clear all
close all
rng(100000)
%% ***************************************** 1. ��ز���
ess_soc=[0 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 1];  % (--)
% Data (Ah): 
ess_max_ah_cap=5;
ess_init_soc=0.9;
ess_r_dis=[0.004178 0.004178 0.004178 0.004178 0.003769 0.003517 0.003296 0.003222 0.003173 0.003197 0.0033 0.003483 0.003732 0.00398 0.00398]*3/5*32; 
ess_r_chg=[0.001393 0.001393 0.001393 0.001393 0.001256 0.001172 0.001099 0.001074 0.001058 0.001066 0.0011 0.001161 0.001244 0.001327 0.001327]*3/5*32; 
ess_voc=[3.124 3.36 3.42 3.479 3.501 3.524 3.547 3.569 3.593 3.616 3.65 3.684 3.72 3.756 3.9]*3*32;
%% **************************************   2.  motor  ����
load('motor.mat')
mc_map_trq=motorMAP(1,2:end);
mc_map_spd=[0;1000;2000;3000;3500;4000;5000;6000;6500;7000;8000]*(2*pi)/60;% (rad/s), speed range of the motor
mc_eff_map =[0.43,0.49,0.58,0.66,0.73,0.70,0.73,0.66,0.58,0.49,0.43;
    0.69,0.73,0.80,0.8440,0.8850,0.70,0.8850,0.8440,0.80,0.73,0.69;
    0.8020,0.83,0.8650,0.8870,0.92,0.70,0.92,0.8870,0.8650,0.83,0.8020;
    0.85,0.8620,0.8820,0.9010,0.9030,0.70,0.9030,0.9010,0.8820,0.8620,0.85;
    0.84,0.8680,0.8820,0.9010,0.9030,0.70,0.9030,0.9010,0.8820,0.8680,0.84;
    0.85,0.87,0.8850,0.8970,0.90,0.70,0.90,0.8970,0.8850,0.87,0.85;
    0.88,0.88,0.8820,0.8940,0.8850,0.70,0.8850,0.8940,0.8820,0.88,0.88;
    0.88,0.88,0.89,0.89,0.8810,0.70,0.8810,0.89,0.89,0.88,0.88;
    0.88,0.88,0.88,0.89,0.88,0.70,0.88,0.89,0.88,0.88,0.88;
    0.88,0.88,0.89,0.89,0.8780,0.70,0.8780,0.89,0.89,0.88,0.88;
    0.88,0.88,0.89,0.89,0.87,0.70,0.87,0.89,0.89,0.88,0.88];
mc_max_trq=[142.2;148.6;146.5;92.4;83.9;71.1;56.9;48.4;43.4;37.7;35.6];  % (N*m)

for i=1:11
    for j=1:11
        if mc_eff_map(i,j)<0.8
            mc_eff_map(i,j)=mc_eff_map(i,j)*0.8;
        else
            mc_eff_map(i,j)=mc_eff_map(i,j)*1.02;
            if mc_eff_map(i,j)>0.99
                mc_eff_map(i,j)=0.99;
            end
        end
    end
    
end
mc_max_gen_trq=-1*mc_max_trq;  % (N*m), estimate
[Tm,wm]=meshgrid(mc_map_trq, mc_map_spd);

% figure
% subplot(211),surf(wm,Tm,mc_eff_map);shading interp;colormap('jet');xlabel ת��rad/s;ylabel ת��N.m;colorbar
% subplot(212),contourf(wm,Tm,mc_eff_map,'ShowText','on');hold on;xlabel ת��rad/s;ylabel ת��N.m;colorbar
% plot(mc_map_spd,mc_max_trq,'k','LineWidth',2);hold on;
% plot(mc_map_spd,mc_max_gen_trq,'k','LineWidth',2);hold on;
%% ****************************************** �Ż��㷨��ʼ
NPar = 2; % �Ż������ĸ���
VarLow =[0 -2]; % ���ٶ�����
VarHigh =[2 0]; % ���ٶ����ޡ�
%% ********************************����
v_base=30;  %Ѳ������
% ax_1=0.2;  %�������ٶ�
% ax_2=-1;%���м��ٶ�
ess_init_soc=0.9;


PopSize =20;  %��������
MaxIterations = 10;%��������
KeepPercent = 20/100;  %���򲻱����
CrossPercent = 40/100; %����������
MutatPercent = 1 - KeepPercent - CrossPercent;
SelectionMode = 1;
KeepNum = round(KeepPercent * PopSize);
CrossNum = round(CrossPercent * PopSize);
if mod(CrossNum,2) ~= 0
    CrossNum = CrossNum + 1;
end
MutatNum = PopSize - KeepNum - CrossNum;
%% ******************��ʼ����Ⱥ
AA=rand(PopSize,NPar);
AB=VarHigh - VarLow;
ABC=[];
for i=1:NPar
    C=AA(:,i)*AB(:,i)+VarLow(1,i);
    ABC=[ABC C];
end
Pop = ABC;
clear AA AB ABC C;
%% ******************��ʼ����Ⱥ��Ӧ�ȼ���
for i=1:PopSize
    ax_1=Pop(i,1);
    ax_2=Pop(i,2);
    sim('ACC_CC_model_best_2017bV2.mdl' );
    error(i)= soc(1)-soc(end);  %���ĵĵ���Խ��Խ��
end
[Cost Indx] = sort(error);
Pop = Pop(Indx,:);
%% ****************************************��ѭ��
MinMat =min(Cost);
for Iter = 1:MaxIterations
    %% **********************************Select Keep
    Pop(1:KeepNum,:) = PSO_Fcn(Pop(1:KeepNum,:));
    %% **********************************����
    SlectedIndexes = SelectParents_Fcn(Cost,CrossNum,SelectionMode);
    CrossPop = [];
    for ii = 1:2:CrossNum
        Par1Index = SlectedIndexes(ii);
        Par2Index = SlectedIndexes(ii+1);
        Par1 = Pop(Par1Index,:);
        Par2 = Pop(Par2Index,:);
        [Off1 , Off2] = CrossOver_fcn(Par1,Par2);
        CrossPop = [CrossPop ; Off1 ; Off2];
    end
    Pop(KeepNum+1:KeepNum+CrossNum,:) = CrossPop;
    %% **********************************����ͻ��
    AA=rand(MutatNum,NPar);
    AB=VarHigh - VarLow;
    ABC=[];
    for i=1:NPar
        C=AA(:,i)*AB(:,i)+VarLow(1,i);
        ABC=[ABC C];
    end
    Pop(KeepNum+CrossNum+1 : end , :)= ABC;
    clear AA AB ABC;
    [m,n] = size(Pop);
    %% **********************************�µ���Ⱥ
    for i=1:n
        jj=1;
        while jj<=m
            if  Pop(jj,i)>VarHigh(1,i)
                Pop(jj,i)=VarHigh(1,i);
            end
            if  Pop(jj,i)<VarLow(1,i)
                Pop(jj,i)=VarLow(1,i);
            end
            jj=jj+1;
        end
    end
    %%****************************************%��Ӧ�ȼ���start
    for i=1:PopSize
    ax_1=Pop(i,1);
    ax_2=Pop(i,2);
    sim('ACC_CC_model_best_2017bV2.mdl' );
    error(i)= soc(1)-soc(end);  %���ĵĵ���Խ��Խ��
    end
    %%****************************************%��Ӧ�ȼ���end
    [Cost Indx] = sort(error);
    Pop = Pop(Indx,:);
    MinMat = [MinMat min(Cost)];
end
%% ***************************************
BestSolution = Pop(1,:)  %���ű���
BestCost = Cost(1)   %������Ӧ��ֵ
plot(MinMat);xlabel ��������;ylabel ��Ӧ��ֵ









