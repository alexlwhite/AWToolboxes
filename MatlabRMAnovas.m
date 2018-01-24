%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeated-measures ANOVAs in MATLAB
% by Alex White 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Requires a matrix of data, ds, in which the last dimension indexes subjects
%% Make some fake data with 10 subjects: 

% 1-way anova: 
% expMeans = [1; 4];
% ds = repmat(expMeans,[1 10]);

% 2-way anova:
% expMeans = [5 5.1; 1 10];
% ds = repmat(expMeans, [1 1 10]);

%Add some labels for the conditions (only for 2-way anova) 
factnames = {'cueCondition','stimulusOrder'};

% 3-way anova 
expMeans = rand(2,2,2)*10; 
expMeans(:,:,2) = expMeans(:,:,2)+5;
ds = repmat(expMeans,[1 1 1 10]); 

%add noise to data: 
ds = ds+randn(size(ds));

%number of subjects: 
nsubj = size(ds,ndims(ds));

%% Specifiy a file to print results to, statsF 
%Something like: 
%statsF=fopen('ANOVARes.txt', 'w');

%But for now, just print to command line 
statsF = 1; 

%% Initialize variables  

%how many experimental variables there are, determines which type of ANOVA
fni = ndims(ds)-1;

g1=zeros(1,numel(ds));
g2=zeros(1,numel(ds));
g3=zeros(1,numel(ds));
g4=zeros(1,numel(ds));
y =zeros(1,numel(ds));

%% ONE-WAY ANOVA 
if fni == 1
    [anvp, stats] = anova_rm(ds',0);
    for sxi=1:size(stats,1)
        for syi=1:size(stats,2)
            if ischar(stats{sxi,syi})
                fprintf(statsF,sprintf('%s\t',stats{sxi,syi}));
            else
                fprintf(statsF,sprintf('%.6f\t',stats{sxi,syi}));
            end
        end
        fprintf(statsF,sprintf('\n'));
    end
    %Note: this function always calls the experimental variable "time"

%% TWO-WAY ANOVA
elseif fni == 2
    dpnum=0;
    for v1=1:size(ds,1)
        for v2=1:size(ds,2)
            for v3=1:size(ds,3)
                dpnum=dpnum+1;
                g1(dpnum)=v1;
                g2(dpnum)=v2;
                g3(dpnum)=v3;
                
                y(dpnum)=ds(v1,v2,v3);
            end
        end
    end
    
    stats=rm_anova2(y,g3,g1,g2,factnames);
    
    
    for sxi=1:size(stats,1)
        for syi=1:size(stats,2)
            if ischar(stats{sxi,syi})
                fprintf(statsF,sprintf('%s\t',stats{sxi,syi}));
            else
                fprintf(statsF,sprintf('%.6f\t',stats{sxi,syi}));
            end
        end
        fprintf(statsF,sprintf('\n'));
    end
    
%% 3-WAY ANOVA
elseif fni == 3
    
    dpnum=0; gnum=0;
    for v1=1:size(ds,1)
        
        gnum=gnum+1;
        for v2=1:size(ds,2)
            for v3=1:size(ds,3)
                for vsubj=1:nsubj
                    dpnum=dpnum+1;
                    y(dpnum)=ds(v1,v2,v3,vsubj);
                    g1(dpnum)=v1;
                    g2(dpnum)=v2;
                    g3(dpnum)=v3;
                    g4(dpnum)=vsubj;
                end
            end
        end
    end
    dataMat = [y' g1' g2' g3' g4'];
    RMAOV33Save(dataMat,0.05,statsF);
end