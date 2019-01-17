%% make boxplot of RO for all basins and each scenario / 30 year period     
    
clear all
plotRO = 1;
plotROfrac = 0;
plotPR = 1;
plotTM = 1;
savePlots = 1;

load 'data/simResults/MoroccoAll22GCMs_15Basins_and16Dams_totals_monthYR.mat'
    
gcms = length(GCMs);

% make monthly mean (remove years)
RObasmn = squeeze(nanmean(RObas,5));
ROdammn = squeeze(nanmean(ROdam,5));
ROdammnOBS = squeeze(nanmean(ROdamOBS,2));
RObasmnOBS = squeeze(nanmean(RObasOBS,2));

% % replace zeros with nan
% Zeros = find(RObasmn == 0);
% RObasmn(Zeros) = nan(size(Zeros));

% make into percent difference from OBS
RObasmnOBSshift = shiftdim(RObasmnOBS,-3);
ROdammnOBSshift = shiftdim(ROdammnOBS,-3);
RObasmnDEL = 1 + (RObasmn - repmat(RObasmnOBSshift,[gcms 3 10 1 1])) ./ repmat(RObasmnOBSshift,[gcms 3 10 1 1]); 
ROdammnDEL = 1 + (ROdammn - repmat(ROdammnOBSshift,[gcms 3 10 1 1])) ./ repmat(ROdammnOBSshift,[gcms 3 10 1 1]);

RO30yr = nan(gcms,3,3,15,12);
for GCM = 1:gcms
for sres = 1:3
    for bas = 1:15
        for yr30 = 1:3
            ROi = squeeze(nanmean(RObasmn(GCM,sres,(yr30-1)*3+2:yr30*3+1,bas,:),3));
            RO30yr(GCM,sres,yr30,bas,:) = ROi;
        end
    end
end
end

RObasOBSan = mean(RObasmnOBS,2);
% replace zeros with nan
Zeros = find(RO30yr == 0);
RO30yr(Zeros) = nan(size(Zeros));
basSTRshort = [

    'Sebou  '
    'Massa  '
    'Draa   '
    'Souss  '
    'Tensift'
    'OumErRb'
    'Loukous'
    'ZizRher'
    'Moulouy'
    'CotMed '
    'CotGueT'
    'Tamri  '
    'CotEssa'
    'CotAtla'
    'Guir   '];

decSTR = ['2030';'2050';'2080'];
sceSTR = ['A1 ';'A1B';'B1 '];
RO30yran = nanmean(RO30yr,5);

f = 0;
LBL = cell(9,1);
LBL{1} = 'A2';
LBL{2} = 'A1B';
LBL{3} = 'B1';
LBL{4} = ' A2 ';
LBL{5} = ' A1B ';
LBL{6} = ' B1 ';
LBL{7} = '  A2  ';
LBL{8} = '  A1B  ';
LBL{9} = '  B1  ';
RO30yr2an = nan(gcms,9,15);
for i = 1:9
    dec3 = ceil(i/3);
    sres = rem(i-1,3)+1;
    RO30yr2an(:,i,:) = squeeze(RO30yran(:,sres,dec3,:));
end

conv = 60*60*24*30.5/1000000;
if plotRO
for bas = 1:15
        BPmat = squeeze(RO30yr2an(:,:,bas))*conv;
        f = f + 1;
        hf = figure(f);
        boxplot(BPmat,LBL);
        ht = title([basSTR{bas},' Runoff   (Dashed Black Line = 1979 - 1990 Mean)']);
        set(ht,'Interpreter','none')
        ylabel('Average Runoff (Millions of Cubic Meters Per Month)')
%         xlabel('Scenario')
        axis([0 10 0 max(nanmax(BPmat))*(1.05)])
        line([0,10],[RObasOBSan(bas)*conv,RObasOBSan(bas)*conv],'linestyle','--','color','k');
                hold on
                text(2,max(max(BPmat))*-0.1,'2030','HorizontalAlignment','center')
                text(2,max(max(BPmat))*-0.063,'-----------------------------','HorizontalAlignment','center');
                text(5,max(max(BPmat))*-0.1,'2050','HorizontalAlignment','center')
                text(5,max(max(BPmat))*-0.063,'-----------------------------','HorizontalAlignment','center');
                text(8,max(max(BPmat))*-0.1,'2080','HorizontalAlignment','center')
                text(8,max(max(BPmat))*-0.063,'-----------------------------','HorizontalAlignment','center');
                hold off
        if savePlots
        saveas(hf,['plots/BoxPlot_Basin_',basSTR{bas},'_Runoff.jpg'])
        end
end
end

RO30yrDEL = nan(gcms,3,3,15,12);
for GCM = 1:gcms
for sres = 1:3
    for bas = 1:15
        for yr30 = 1:3
            ROi = squeeze(nanmean(RObasmnDEL(GCM,sres,(yr30-1)*3+2:yr30*3+1,bas,:),3));
            RO30yr(GCM,sres,yr30,bas,:) = ROi;
        end
    end
end
end

RO30yranDEL = nanmean(RO30yrDEL,5);
% % replace zeros with nan
% Zeros = find(RO30yranDEL == 0);
% RO30yranDEL(Zeros) = nan(size(Zeros));

RO30yr2anDEL = nan(gcms,9,15);
for i = 1:9
    dec3 = ceil(i/3);
    sres = rem(i-1,3)+1;
    RO30yr2anDEL(:,i,:) = squeeze(RO30yranDEL(:,sres,dec3,:));
end

% plot boxpolots of % difference
if plotROfrac
for bas = 1:15
        BPmat = squeeze(RO30yr2anDEL(:,:,bas)) .* 100;
        f = f + 1;
        hf = figure(f);
        boxplot(BPmat,LBL);
        ht = title([basSTR{bas},' Change in Runoff   (Dashed Black Line = 1979 - 1990 Mean)']);
        set(ht,'Interpreter','none')
        ylabel('Runoff (% Difference from Observed)')
%         xlabel('Scenario')

        rng = zeros(3,1);
        rng(3) = max(nanmax(BPmat)) - min(nanmin(BPmat));
        rng(1) = min((min(nanmin(BPmat)) - rng(3)*0.05),rng(3)*-0.1);
        rng(2) = max(nanmax(BPmat)) + rng(3)*0.05;
        axis([0 10 rng(1) rng(2)])
        
        line([0,10],[0,0],'linestyle','--','color','k');
                hold on
                text(2,rng(1) - rng(3)*0.1,'2030','HorizontalAlignment','center')
                text(2,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                text(5,rng(1) - rng(3)*0.1,'2050','HorizontalAlignment','center')
                text(5,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                text(8,rng(1) - rng(3)*0.1,'2080','HorizontalAlignment','center')
                text(8,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                hold off
        if savePlots
        saveas(hf,['plots/BoxPlot_Basin_',basSTR{bas},'_Runoff_perchange.jpg'])
        end
end
end
    
%% plot Precip boxplots
clear clear RObas RO30yran RO30yr RObasmn
load 'data/simResults/MoroccoAll22GCMs_15Basins_and16Dams_meanPRCPandTEMP_monthYR.mat'

gcms = length(GCMs);

PR30yr = nan(gcms,3,3,15,12);
for GCM = 1:gcms
for sres = 1:3
    for bas = 1:15
        for yr30 = 1:3
            PRi = squeeze(nanmean(PRbas(GCM,sres,(yr30-1)*3+2:yr30*3+1,bas,:),3));
            PR30yr(GCM,sres,yr30,bas,:) = PRi;
        end
    end
end
end

% replace zeros with nan
Zeros = find(PR30yr == 0);
PR30yr(Zeros) = nan(size(Zeros));

PR30yran = mean(PR30yr,5);

PR30yr2an = nan(gcms,9,15);
for i = 1:9
    dec3 = ceil(i/3);
    sres = rem(i-1,3)+1;
    PR30yr2an(:,i,:) = squeeze(PR30yran(:,sres,dec3,:));
end

% plot boxpolots of % difference in Precip
if plotPR
for bas = 1:15
        BPmat = (squeeze(PR30yr2an(:,:,bas)) - 1) * 100;
        f = f + 1;
        hf = figure(f);
        boxplot(BPmat,LBL);
        ht = title([basSTR{bas},' Change in Precipitation   (Dashed Black Line = 1979 - 1990 Mean)']);
        set(ht,'Interpreter','none')
        ylabel('Precipitation (% Difference from Observed)')
%         xlabel('Scenario')
        rng = zeros(3,1);
        rng(3) = max(nanmax(BPmat)) - min(nanmin(BPmat));
        rng(1) = min((min(nanmin(BPmat)) - rng(3)*0.05),rng(3)*-0.1);
        rng(2) = max(nanmax(BPmat)) + rng(3)*0.05;
        axis([0 10 rng(1) rng(2)])
        
        line([0,10],[0,0],'linestyle','--','color','k');
                hold on
                text(2,rng(1) - rng(3)*0.1,'2030','HorizontalAlignment','center')
                text(2,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                text(5,rng(1) - rng(3)*0.1,'2050','HorizontalAlignment','center')
                text(5,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                text(8,rng(1) - rng(3)*0.1,'2080','HorizontalAlignment','center')
                text(8,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                hold off
        if savePlots
        saveas(hf,['plots/BoxPlot_Basin_',basSTR{bas},'_Precip_perchange.jpg'])
        end
end
end


%% plot TM boxplots

load 'data/simResults/MoroccoAll22GCMs_15Basins_and16Dams_meanPRCPandTEMP_monthYR.mat'

gcms = length(GCMs);

TM30yr = nan(gcms,3,3,15,12);
for GCM = 1:gcms
for sres = 1:3
    for bas = 1:15
        for yr30 = 1:3
            TMi = squeeze(nanmean(TMbas(GCM,sres,(yr30-1)*3+2:yr30*3+1,bas,:),3));
            TM30yr(GCM,sres,yr30,bas,:) = TMi;
        end
    end
end
end

% replace zeros with nan
Zeros = find(TM30yr == 0);
TM30yr(Zeros) = nan(size(Zeros));
    
TM30yran = mean(TM30yr,5);

TM30yr2an = nan(gcms,9,15);
for i = 1:9
    dec3 = ceil(i/3);
    sres = rem(i-1,3)+1;
    TM30yr2an(:,i,:) = squeeze(TM30yran(:,sres,dec3,:));
end

% plot boxpolots of % difference in Precip
if plotTM
for bas = 1:15
        BPmat = (squeeze(TM30yr2an(:,:,bas)));
        f = f + 1;
        hf = figure(f);
        boxplot(BPmat,LBL);
        ht = title([basSTR{bas},' Change in Temperature   (Dashed Black Line = 1979 - 1990 Mean)']);
        set(ht,'Interpreter','none')
        ylabel('Temperature (Difference from Observed (deg C))')
%         xlabel('Scenario')
        rng = zeros(3,1);
        rng(3) = max(nanmax(BPmat)) - min(nanmin(BPmat));
        rng(1) = min((min(nanmin(BPmat)) - rng(3)*0.05),rng(3)*-0.1);
        rng(2) = max(nanmax(BPmat)) + rng(3)*0.05;
        axis([0 10 rng(1) rng(2)])
        
        line([0,10],[0,0],'linestyle','--','color','k');
                hold on
                text(2,rng(1) - rng(3)*0.1,'2030','HorizontalAlignment','center')
                text(2,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                text(5,rng(1) - rng(3)*0.1,'2050','HorizontalAlignment','center')
                text(5,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                text(8,rng(1) - rng(3)*0.1,'2080','HorizontalAlignment','center')
                text(8,rng(1) - rng(3)*0.063,'-----------------------------','HorizontalAlignment','center');
                hold off
        if savePlots
        saveas(hf,['plots/BoxPlot_Basin_',basSTR{bas},'_Temp_perchange.jpg'])
        end
end
end   
    
    
    
    