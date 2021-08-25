function ExcitationContraction
%major overhaul of segmenting data on V16, V15 and previous used peak
%detection, rejection of peaks and other ways to detect contractions. this
%is cleaner and should have been done from the start. -AS 8/6/18

cd 'C:\Users\79829\Documents\Data\OrmerodScibelli\';
[file,path] = uigetfile('*.mat');
cd(path);
load(file);

sf = 2000;       %sample freq., prior to Aug2 sf=2000, post Aug3 sf=5000.
%%Raw data filtering 
%simple average
simpleavg = movmean(data(1,:),10);   %Smoothing filter design, normally 100
simpleavg2 = movmean(data(1,:),1000);   %Smoothing filter design

%% Trial segmenting
%
thresh = 1;                   
dataend = (length(data)-(sf*5));
positrig = find(data(2,(sf*1):dataend) >= thresh);        %if stim at start, change to (sf*1):dataend, normally 1:dataend
%positrig(1,:) = positrig(1,:)+(sf*5),      %add this back in if there is a
%stim very close to the beginning of the trial.
%find(data(2,((sf*5):dataend)) etc.

intraburst = [1];
intraburst = [intraburst,(find(diff(positrig) >= (sf*1.2)))];
intraburst = [intraburst, length(positrig)];
intraburst = unique(intraburst);
        
Trialslimit(1,1) = positrig(intraburst(1));
if (positrig(intraburst(2)))-positrig(intraburst(1)) >= (sf*1.2)
    Trialslimit(1,2) = positrig(intraburst(2));
else
    Trialslimit(1,2) = positrig(intraburst(1));
end

i=1;
for i = 2: length(intraburst)-1
    Trialslimit(i,1) = (positrig(intraburst(i)+1));
    if ((positrig(intraburst(i+1)))-(positrig(intraburst(i)+1))) >= (sf*1.2)    %what is this doing? AS7/2/19
        Trialslimit(i,2) = (positrig(intraburst(i+1)));
    else
        Trialslimit(i,2) = positrig(intraburst(i)+1);
    end
end

i = 1;
for i = 1:size(Trialslimit)
    Trials{i,1} = simpleavg(Trialslimit(i,1)-(sf*0.5):Trialslimit(i,2)+(sf*3));  
    Trials{i,2} = data(2,(Trialslimit(i,1)-(sf*0.5):Trialslimit(i,2)+(sf*3)));  
    baseline = mean(Trials{i,1}(sf*0.01:sf*0.1));
    Trials{i,1} = (Trials{i,1}-(baseline));
    [Trials{i,3},peakcross] = max(Trials{i,1});          %find the max voltage
    Trials{i,3} = round(Trials{i,3},3);            %round to nearest thousandth
    stimtrace = find(Trials{i,2}(1,1:length(Trials{i,1}))>thresh);
    taucross1 = find(Trials{i,1}(1,1:length(Trials{i,1}))>Trials{i,3}*(1-(1/exp(1))));                %tau rise
    taucross2 = find(Trials{i,1}(1,1:length(Trials{i,1}))>Trials{i,3}*((1/exp(1))));                %tau decay
    stimrange = (stimtrace(end)+(0.1*sf));
    if peakcross <= stimrange
        Trials{i,4} = ((taucross1(1)-stimtrace(1))/(sf*0.001));          %divided at the end so that it would be in milliseconds.
        Trials{i,5} =  ((taucross2(end)-stimtrace(end))/(sf*0.001));                  %tau decay
    else
        Trials{i,4} = ((taucross1(1)-stimtrace(1))/(sf*0.001));          %divided at the end so that it would be in milliseconds.
        Trials{i,5} =  ((taucross2(end)-peakcross)/(sf*0.001));                  %tau decay
    end
end


%% Grouping
%plot data trace, with left click at points between stim changes, take 0 to
%first click, 1st to second, etc. n-1 to end. enter to end 
groupquery = input('\nDo you want to group data?(y or n)\n**Push ENTER when done**\n','s');
if groupquery == 'y'
    figure; plot(data(1,:)); hold on;           %to zoom in and acurately mark the low and long thresholds, change to (1,1:800000).
    figtitle = ['Mark between groups. Press Enter to finish.'];
    title(figtitle);
    grouppts = getpts;
    close
    i = 1;
    tempgroup = find(Trialslimit(:,1) <= grouppts(1));
    for i = 1:length(tempgroup)
        groups{1,i}(1,:) = Trials{tempgroup(i),1};
        groups{1,i}(2,:) = Trials{tempgroup(i),2};
        groupsdata{1,1}(i,1) = Trials{tempgroup(i),3};
        groupsdata{1,1}(i,2) = Trials{tempgroup(i),4};
        groupsdata{1,1}(i,3) = Trials{tempgroup(i),5};
    end
elseif groupquery(1) == 'n'    
end

if groupquery(1) == 'y'
    i = 1;
    for i = 1:length(grouppts)-1
        tempgroup = find(Trialslimit(:,1) >= grouppts(i) & Trialslimit(:,1) <= grouppts(i+1));
        ii = 1;
        for ii = 1:length(tempgroup)
            groups{i+1,ii}(1,:) = Trials{tempgroup(ii),1};
            groups{i+1,ii}(2,:) = Trials{tempgroup(ii),2};
            groupsdata{i+1,1}(ii,1) = Trials{tempgroup(ii),3};
            groupsdata{i+1,1}(ii,2) = Trials{tempgroup(ii),4};
            groupsdata{i+1,1}(ii,3) = Trials{tempgroup(ii),5};
        end
    end
elseif groupquery(1) == 'n'
end

if groupquery(1) == 'y'
    tempgroup = find(Trialslimit(:,1) >= grouppts(i+1));
    i = i+1;
    for ii = 1:length(tempgroup)
        groups{i+1,ii}(1,:) = Trials{tempgroup(ii),1};
        groups{i+1,ii}(2,:) = Trials{tempgroup(ii),2};
        groupsdata{i+1,1}(ii,1) = Trials{tempgroup(ii),3};
        groupsdata{i+1,1}(ii,2) = Trials{tempgroup(ii),4};
        groupsdata{i+1,1}(ii,3) = Trials{tempgroup(ii),5};
    end
elseif groupquery(1) == 'n'
end
   

%% Averaging
if groupquery(1) == 'y'
    i = 1;
    groupsize = size(groups);
    for i = 1:groupsize(1)
        ii = 1;
        for ii = 1:sum(~cellfun('isempty',groups(i,:)))
            groupslength(ii) = length(groups{i,ii});
        end
        mingroupslength = min(groupslength)-(sf*0.05);
        ii = 1;
        for ii = 1:sum(~cellfun('isempty',groups(i,:)))
            groupsmatrix{i,1}(ii,:) = groups{i,ii}(1,1:mingroupslength);
        end
        clear groupslength mingroupslength;
    end  
elseif groupquery(1) == 'n'
end

%% Confidence interval
if groupquery(1) == 'y'
    i = 1;
    groupsmatrixsize = size(groupsmatrix);
    for i = 1:groupsmatrixsize(1)
        tempgroupsmatrixsize = size(groupsmatrix{i,1});
        if tempgroupsmatrixsize(1) == 1
        elseif tempgroupsmatrixsize(1) > 30      %determine what to use at z or t. base off cycle number.
            z = 1.96;             %use z table which is static i guess.
            t = 0;
        else
            z = 0;
            ttable = [12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042];
            t = ttable(1,(tempgroupsmatrixsize(1))-1);      %pulled from t-table look up for confidence interval with under 25 = n.
        end
        if tempgroupsmatrixsize(1) == 1
        else
            groupsMEAN{i} = mean(groupsmatrix{i,1});
            groupsSTD{i} = std(groupsmatrix{i,1});
            groupsCI{1,i} = groupsMEAN{1,i}+(z+t).*(groupsSTD{1,i}/sqrt(tempgroupsmatrixsize(1)));
            groupsCI{2,i} = groupsMEAN{1,i}-(z+t).*(groupsSTD{1,i}/sqrt(tempgroupsmatrixsize(1)));
        end    
    end
elseif groupquery(1) == 'n'
end

%groupsMEANratio=groupsMEAN;
groupsdatasize=size(groupsdata);
i=1;
for i=1:(groupsdatasize(1,1))
    tmpmax = mean(groupsdata{i,1});
    %math divide all values by max force, times by 100. get 0-100 scale of
    %data.
    groupsMEANratio{i,1}=groupsMEAN{1,i} ./ tmpmax(1,1);
    %clear tmpmax
end

%% Plotting
% for plotting all trials in order before sorting or grouping
if groupquery(1) == 'y'
    i = 1;
    ii = 1;
    %figure, hold on;           %figure overlay, comment out to plot groups independant, comment out from for loop
    for i = 1:groupsize(1)
        if size(groupsmatrix{i,1}) > 1
            tempgroupsmatrixsize = size(groupsmatrix{i,1});
            figure; hold on;                %substitute this for figure hold before for loop
            figtitle = ['Group #' num2str(i)];
            title(figtitle);
            % Title the group figure
            yy1 = groupsCI{1,i};
            yy2 = groupsCI{2,i};
            x = linspace(1,length(yy1),length(yy1));
            fillcolor = [.75 .75 .75];
            fill([x fliplr(x)],[yy1 fliplr(yy2)],fillcolor);        %,'edgecolor','none' after fillcolor
            plot(groupsMEAN{i},'k','LineWidth',1)
            for ii = 1:tempgroupsmatrixsize(1)     %also comment this out
            %to remove raw data 
                plot(groupsmatrix{i,1}(ii,:));
            end  
        else
            figure; hold on;
            figtitle = ['Group #' num2str(i)];
            title(figtitle);
            plot(groupsmatrix{i,1})
        end
    end
elseif groupquery(1) == 'n'
    figure; hold on;
    i=1;
    for i =1:(length(Trials))
        plot(Trials{i,1});
        pause(0.25);
    end
end

%Normalized average plotting
figure, hold on;
i=1;
for i = 1:(groupsdatasize(1,1))
    plot(groupsMEANratio{i,1});
end

%% Saving routine
if groupquery(1) == 'y'
    uisave({'data','sf','Trials','groupsdata','groupsMEANratio','Trialslimit','groups','groupsmatrix','groupsMEAN','groupsCI','file','path'}, '_processed')
elseif groupquery(1) == 'n'
    uisave({'data','sf','Trials','Trialslimit','file','path'}, '_processed')
end
close all

%home,disp('Press enter to end'), pause;
