
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%                            Introduction                                 %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Hey there, 
%
% This script performs all analyses carried out in the Neuroimage Clinical
% paper 'Quasi-periodic patterns of brain activity in individuals with
% attention-deficit/hyperactivity disorder' from start to finish. If you
% have any questions about the methods, please feel free to contact me via
% email. I encourage you to publish the code for your methods alongisde
% your papers as well. 
%
% Cheers,
% Anzar Abbas, PhD
% anzar@anzarabbas.com

for setting_dirs = 1:1
    if ismac == 1
     dir_group1 = '/Volumes/keilholz-lab/Anzar/11-paper_adhd/Data_Control';
     dir_group2 = '/Volumes/keilholz-lab/Anzar/11-paper_adhd/Data_ADHD';
     dir_results = '/Volumes/keilholz-lab/Anzar/11-paper_adhd/Results';
     dir_atlas = '/Volumes/keilholz-lab/Anzar/11-paper_adhd/Atlas';
     dir_groups = {dir_group1,dir_group2};
    else
     dir_group1 = '/keilholz-lab/Anzar/11-paper_adhd/Data_Control';
     dir_group2 = '/keilholz-lab/Anzar/11-paper_adhd/Data_ADHD';
     dir_results = '/keilholz-lab/Anzar/11-paper_adhd/Results';
     dir_atlas = '/keilholz-lab/Anzar/11-paper_adhd/Atlas';
     dir_groups = {dir_group1,dir_group2};
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%                   Data Management & Preprocessing                       %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% I am going to download data from the ADHD 200 dataset, which is part of
% the 1000 Functional Connectome Project. Since I will have to use data
% from multiple sites, it would be preferable if the scan parameters are
% similar for the data collected across sites. Hence, I will be downloading
% data from the New York University, Peking University, and Neuroimage
% datasets. At the very least, the TR for all three of those datasets was 
% 2 seconds.

for downloading_data = 1:1
    
    adhd200 = ' ftp://www.nitrc.org/fcon_1000/htdocs/indi/adhd200/sites/';
    % This is the website where the ADHD 200 datasets are stored
    
    cd /keilholz-lab/Anzar/11-paper_adhd/RawData/NewYorkUniversity
    system(['wget',adhd200,'/nyu/data/NYU.001.002.part1.tar.gz']);
    system(['wget',adhd200,'/nyu/data/NYU.001.002.part2.tar.gz']);
    system(['wget',adhd200,'/nyu/data/NYU.001.001.part3.tar.gz']);
    system(['wget',adhd200,'/nyu/data/NYU.001.001.part4.tar.gz']);
    system('tar xfv NYU.001.002.part*')
    % Downloading and unpacking the NYU dataset

    cd /keilholz-lab/Anzar/11-paper_adhd/RawData/PekingUniversity
    system(['wget',adhd200,'/beijing/data/Peking_1.001.003.tar.gz']);
    system(['wget',adhd200,'/beijing/data/Peking_2.001.001.tar.gz']);
    system(['wget',adhd200,'/beijing/data/Peking_3.001.001.tar.gz']);
    system('tar xfv Peking_3.001.00*');
    % Downloading the Peking University dataset
    
    cd /keilholz-lab/Anzar/11-paper_adhd/RawData/NeuroImage
    system(['wget',adhd200,'nijmegen/data/NeuroIMAGE.001.001.tar.gz'])
    system(['wget',adhd200,'nijmegen/data/NeuroIMAGE_TestRelease.001.002.tar.gz'])
    system('tar xfv NeuroIMAGE_TestRelease.001.002.tar.gz');
    system('tar xfv NeuroIMAGE.001.001.tar.gz');
    % Downloading the  NeuroImage dataset
        
end

% ----------------------------------------------------------------------- %

% Next, I need to organize the data I just downloaded. Each dataset came
% with an excel file that has phenotypic data for the individuals scanned.
% I want to separate the controls from the individuals with ADHD. I am also
% only keeping those with 'combined' ADHD. The diagnosis of all individuals
% is saved in the 'dx' column in the excel sheets and I have manually
% extracted that column. Additionally, I am only keeping the first
% resting-state scan from each individual and that too if it passed the
% data quality check from the collectors of the data, saved as 'qc.'

for organizing_data = 1:1
    
    cd /keilholz-lab/Anzar/11-paper_adhd/RawData/NewYorkUniversity
    load('dx.mat'), load('qc.mat'), qc = qc(:,1);
    % I manually extracted this information from the phenotype excel sheets
    % that came with the dataset. dx.mat has the diagnostic information for
    % each individual and qc.mat has the quality check information for the
    % first resting-state scan of each individual.
    cd NYU
    files = dir(pwd);
    subjects = {files.name}';
    directories = [files.isdir]';
    directories(1:2) = 0;
    subjects = subjects(directories);
    for i = 1:length(subjects)
        if dx(i) == 0 && qc(i) == 1
            system(['cp -r ',subjects_control{i},' /keilholz-lab/', ...
                'Anzar/11-paper_adhd/Data_Control/subject_nyu', ...
                subjects_control{i}]);
        end
        if dx(i) == 1 && qc(i) == 1
            system(['cp -r ',subjects_adhd{i},' /keilholz-lab/Anzar/', ...
                '11-paper_adhd/Data_ADHD/subject_nyu',subjects_adhd{i}]);
        end
    end
    % Copying all the New York University control and ADHD scans from the
    % raw data folder to the folders where all the analysis will eventually
    % be carried out

    cd /keilholz-lab/Anzar/11-paper_adhd/RawData/PekingUniversity
    load('dx.mat'), load('qc.mat');
    cd Peking
    files = dir(pwd);
    subjects = {files.name}';
    directories = [files.isdir]';
    directories(1:2) = 0;
    subjects = subjects(directories);
    for i = 1:length(subjects)
        if dx(i) == 0 && qc(i) == 1
            system(['cp -r ',subjects_control{i},' /keilholz-lab/', ...
                'Anzar/11-paper_adhd/Data_Control/subject_peking', ...
                subjects_control{i}]);
        end
        if dx(i) == 1 && qc(i) == 1
            system(['cp -r ',subjects_adhd{i},' /keilholz-lab/Anzar/', ...
                '11-paper_adhd/Data_ADHD/subject_peking', ...
                subjects_adhd{i}]);
        end
    end
    % Copying all the Peking University control and ADHD scans from the raw
    % data folder to the folders where all the analysis will eventually be
    % carried out
    
    cd /keilholz-lab/Anzar/11-paper_adhd/RawData/NeuroImage
    load('dx.mat'), load('qc.mat');
    cd NeuroIMAGE
    files = dir(pwd);
    subjects = {files.name}';
    directories = [files.isdir]';
    directories(1:2) = 0;
    subjects = subjects(directories);
    for i = 1:length(subjects)
        if dx(i) == 0 && qc(i) == 1
            system(['cp -r ',subjects_control{i},' /keilholz-lab/', ...
                'Anzar/11-paper_adhd/Data_Control/subject_neuroimage', ...
                subjects_control{i}]);
        end
        if dx(i) == 1 && qc(i) == 1
            system(['cp -r ',subjects_adhd{i},' /keilholz-lab/Anzar/', ...
                '11-paper_adhd/Data_ADHD/subject_neuroimage', ...
                subjects_adhd{i}]);
        end
    end
    % Copying all the NeuroImage control and ADHD scans from the raw
    % data folder to the folders where all the analysis will eventually be
    % carried out
    
end

% ----------------------------------------------------------------------- %

% Next, I am organizing all the files in a way that is understandable and
% works with the automated preprocessing pipeline that is coming up. This
% is how all my datasets are organized (as is evident by other projects)
% and makes following the steps that were carried out on the data easy to
% follow for anyone that's looking at it.

for cleaning_up_folders = 1:2
    cd (dir_groups{cleaning_up_folders})
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        fprintf([subjects{i},'\n'])
        system('cp session_1/anat_1/mprage* ./t1.nii.gz');
        system('cp session_1/rest_1/rest.nii.gz ./f01.nii.gz');
        system('rm -rf session*')
        cd ..
    end
    % Cleaning up the folders for all datasets
    
    % After doing this, I manually only kept as many control scans from
    % each dataset as there were ADHD scans. This will make sure that the
    % data points from both groups are equal and I'm not wasting time
    % analyzing data points I won't actually end up using in the stats.
    
end

% ----------------------------------------------------------------------- %

% Finally, I am going to run my automated preprocessing pipeline to get the
% data ready for analysis. This preprocessing pipeline is available in the
% same directory as this script and also on the Keilholz lab website. For
% preprocessing, I am first cutting all the functional scan lengths to the
% lowest common denominator (170 timepoints). The anatomical scans will be
% registered to MNI atlas space and tissue segmented. The functional scans
% will be slice time corrected, motion corrected, registered to MNI,
% spatially smoothed, temporally filtered, white matter and CSF signal
% regressed, global signal regressed, and z-scored. FYI, everything beyond
% here was also done without global signal regression. Those results are in
% the supplementary figures in the paper.

for preprocessing = 1:2
    cd (dir_groups{preprocessing})
    atlasdir = '/keilholz-lab/Anzar/Tools/Preprocess/Atlases/Human';
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        system('fslroi f01.nii.gz  f01.nii.gz  0 170')
        anz_preprocess(pwd,atlasdir,2,2,10,[1 2],6,[.01 .08],1,1,1)
        cd ..
    end
    % Preprocessing all of the data the same way it was done in Abbas
    % et al. 2018
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%                 Comparing DMN and TPN across groups                     %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% The Brainnetome atlas was downloaded (Fan et al. 2016) and the posterior
% cingulate cortex was isolated and saved. Here, I am going to get the
% voxels in the brain that are highly correlated anticorrelated with the
% PCC's timecourse and call them the DMN and TPN respectively. Then, I am
% going to find out what the names of those regions are and compare them
% across the control and ADHD group.

for getting_networks = 1:2
    
    cd (dir_groups{getting_networks})
    fprintf(['\n',whatsthetime,' - Getting DMN and TPN for group ', ...
        getting_networks,'\n'])
    
    f = [];
    subjects_to_use = [1:10,24:33,81:90];
    % This combination of subjects ensures that I am getting ten scans from
    % each dataset, so that  one dataset does not bias the results.
    
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects_to_use)
        cd (subjects{subjects_to_use(i)})
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D.mat']);
        f = cat(3,f,func);
        cd ..
    end
    % Concatenating data from 30 subjects. This should be pretty heavy on
    % the memory, and I don't think I can do more than 30 subjects with the
    % memory limitations I have on the machine I am using for this project.
    
    cd resources
    load('roi_pcc.mat')
    [dmn,tpn,corrmap,thresh_dmn,thresh_tpn] = get_ac_networks(pcc,f,10);
    save('networks','dmn','tpn');
    % Getting maps of DMN and TPN by extracting 10% of the most and least
    % PCC-correlated voxels in the brain. Also, thresh_dmn and thresh_tpn
    % are the 90th and 10th percentile values in corrmap that were used as
    % cutoffs to create the DMN and TPN binary masks.
    
    cd (dir_atlas)
    load('roi_atlas.mat'), load('labels.mat')
    num_rois = length(labels);
    % Loading 273 ROI binary masks and their labels from the Brainnetome
    % atlas acquired from Fan et al. (2016)
    
    dmn_rois = cell(1,7);
    tpn_rois = cell(1,7);
    % These arrays will hold information on the contents of the DMN and TPN
    % masks, including the ROI label number it is in the Brainnetome atlas,
    % the number of  voxels included in the networks in each ROI, and the
    % correlation values for each ROI. This information will be compiled
    % into Supplementary Table 1 and will be very helpful in the analysis.
    
    dmn_roi_count = 0;
    tpn_roi_count = 0;
    for i = 1:num_rois
        
        dmn_roi_test = rois{i} .* dmn;
        % This matrix will have twos wherever the DMN overlapped with the 
        % ROI that is being currently looked at
        if length(find(dmn_roi_test)) > 0 %#ok<ISMT>
            dmn_roi_count = dmn_roi_count + 1;
            % Updating the count
            dmn_rois{dmn_roi_count,1} = i;
            % Saving the ROI index 
            dmn_rois{dmn_roi_count,2} = labels{i};
            % Saving the ROI  name
            dmn_rois{dmn_roi_count,3} = length(find(dmn_roi_test));
            % Saving how many voxels of the ROI are in the DMN
            dmn_roi_corrmap = corrmap .* dmn_roi_test;
            dmn_roi_corrmap(dmn_roi_corrmap==0) = NaN;
            dmn_rois{dmn_roi_count,5} = mean(dmn_roi_corrmap(:),'omitnan');
            % Saving the mean correlation for those voxels
        end
        
        tpn_roi_test = rois{i} + tpn;
        % This matrix will have twos wherever the DMN overlapped with the 
        % ROI that is being currently looked at
        if length(find(tpn_roi_test)) > 0 %#ok<ISMT>
            tpn_roi_count = tpn_roi_count + 1;
            % Updating the count
            tpn_rois{tpn_roi_count,1} = i;
            % Saving the ROI index 
            tpn_rois{tpn_roi_count,2} = labels{i};
            % Saving the ROI  name
            tpn_rois{tpn_roi_count,3} = length(find(tpn_roi_test));
            % Saving how many voxels of the ROI are in the DMN
            tpn_roi_corrmap = corrmap .* tpn_roi_test;
            tpn_roi_corrmap(tpn_roi_corrmap==0) = NaN;
            tpn_rois{tpn_roi_count,5} = mean(tpn_roi_corrmap(:),'omitnan');
            % Saving the mean correlation for those voxels
        end
        
    end
    % For the predefined arrays, filling in the ROI index, the ROI name,
    % the number of voxels in the DMN/TPN for that ROI, and the mean
    % correlation with the PCC for those voxels. The rest w ill be filled
    % out in the next few steps.
    
    col2 = dmn_rois(:,2);
    for i = 1:length(col2), col2{i} = col2{i}(1:10); end
    [~,u,~] = unique(col2,'stable');
    u = cat(1,u,length(col2)+1);
    for i = 1:length(u)-1
        dmn_rois{u(i),7} = zeros(size(rois{1}));
        dmn_rois{u(i),4} = 0;
        dmn_rois{u(i),6} = 0;
        for j = u(i):u(i+1)-1
            dmn_rois{u(i),4} = dmn_rois{u(i),4} + dmn_rois{j,3};
            dmn_rois{u(i),6} = dmn_rois{u(i),6} + dmn_rois{j,5};
            voxels_to_add = dmn .* rois{dmn_rois{j,1}};
            dmn_rois{u(i),7} = dmn_rois{u(i),7} + voxels_to_add;
        end
        dmn_rois{u(i),6} = dmn_rois{u(i),6}/length(u(i):u(i+1)-1);
    end
    % Filling in the rest of the cell array. This includes the total number
    % of voxels in each brain region, the mean correlation with the PCC for
    % the brain region, and the how much DMN the brain region includes.
    
    col2 = tpn_rois(:,2);
    for i = 1:length(col2), col2{i} = col2{i}(1:10); end
    [~,u,~] = unique(col2,'stable');
    u = cat(1,u,length(col2)+1);
    for i = 1:length(u)-1
        tpn_rois{u(i),7} = zeros(size(rois{1}));
        tpn_rois{u(i),4} = 0;
        tpn_rois{u(i),6} = 0;
        for j = u(i):u(i+1)-1
            tpn_rois{u(i),4} = tpn_rois{u(i),4} + tpn_rois{j,3};
            tpn_rois{u(i),6} = tpn_rois{u(i),6} + tpn_rois{j,5};
            voxels_to_add = tpn .* rois{tpn_rois{j,1}};
            tpn_rois{u(i),7} = tpn_rois{u(i),7} + voxels_to_add;
        end
        tpn_rois{u(i),6} = tpn_rois{u(i),6}/length(u(i):u(i+1)-1);
    end
    % Filling in the rest of the cell array. This includes the total number
    % of voxels in each brain region, the mean correlation with the PCC for
    % the brain region, and the how much DMN the brain region includes.
    
    cd (dir_results)
    if getting_networks == 1
        save('1-dmn_tpn_group1','dmn','tpn','corrmap','thresh_dmn',  ...
            'thresh_tpn')
        corrmap_dmn_and_tpn = corrmap .* dmn + corrmap .* tpn;
        unslice_save(corrmap_dmn_and_tpn,[91 109 91], ...
            '3-dmn_tpn_group1.nii.gz',[2 2 2],[45 65 37]);
        save('5-dmn_tpn_info_group1','dmn_rois','tpn_rois');
        dmn_rois_g1 = dmn_rois;
        tpn_rois_g1 = tpn_rois;
        % Saving the results for the control group and saving some
        % variables to be saved alongside the ADHD group
    else
        save('2-dmn_tpn_group2','dmn','tpn','corrmap','thresh_dmn',  ...
            'thresh_tpn')
        corrmap_dmn_and_tpn = corrmap .* dmn + corrmap .* tpn;
        unslice_save(corrmap_dmn_and_tpn,[91 109 91], ...
            '4-dmn_tpn_group2.nii.gz',[2 2 2],[45 65 37]);
        save('6-dmn_tpn_info_group2','dmn_rois','tpn_rois');
        dmn_rois_g2 = dmn_rois;
        tpn_rois_g2 = tpn_rois;
        % Saving the results for the ADHD group and saving some variables
        % for further analysis below. What's left to do pretty much is
        % create the ROI atlas I'll be using for the functional
        % connectivity analysis. I haven't explained this yet but I will
        % when I get to the comparing FC matrices part of the paper.
        
        rois_count = 0;
        rois_dmntpn = cell(1,2);
        % This will be the cell that will hold all DMN and TPN ROIs
        
        dmn_rois = cat(1,dmn_rois_g1,dmn_rois_g2);
        dmn_rois = sortrows(dmn_rois,1);
        dmn_lbls = dmn_rois(:,2);
        for i = 1:size(dmn_lbls,1)
            dmn_lbls{i} = dmn_lbls{i}(1:10);
        end
        [~,u,~] = unique(dmn_lbls,'stable');
        u = cat(1,u,length(dmn_rois(:,2))+1);
        for i = 1:length(u)-1
            rois_count = rois_count + 1;
            rois_dmntpn{rois_count,1} = dmn_rois{u(i),2};
            rois_dmntpn{rois_count,2} = zeros(1090,910);
            dmn_roi = dmn_rois(u(i):u(i+1)-1,:);
            for j = 1:size(dmn_roi,1)
               if length(find(dmn_roi{j,7})) > 0
                   rois_dmntpn{rois_count,2} = ...
                       rois_dmntpn{rois_count,2} + dmn_roi{j,7};
               end
            end
            rois_dmntpn{rois_count,2} = ...
                double(logical(rois_dmntpn{rois_count,2}));
        end
        % Adding the DMN ROIs to the atlas
        
        tpn_rois = cat(1,tpn_rois_g1,tpn_rois_g2);
        tpn_rois = sortrows(tpn_rois,1);
        tpn_lbls = tpn_rois(:,2);
        for i = 1:size(tpn_lbls,1)
            tpn_lbls{i} = tpn_lbls{i}(1:10);
        end
        [~,u,~] = unique(tpn_lbls,'stable');
        u = cat(1,u,length(tpn_rois(:,2))+1);
        for i = 1:length(u)-1
            rois_count = rois_count + 1;
            rois_dmntpn{rois_count,1} = tpn_rois{u(i),2};
            rois_dmntpn{rois_count,2} = zeros(1090,910);
            tpn_roi = tpn_rois(u(i):u(i+1)-1,:);
            for j = 1:size(tpn_roi,1)
               if length(find(tpn_roi{j,7})) > 0
                   rois_dmntpn{rois_count,2} = ...
                       rois_dmntpn{rois_count,2} + tpn_roi{j,7};
               end
            end
            rois_dmntpn{rois_count,2} = ...
                double(logical(rois_dmntpn{rois_count,2}));
        end
        % Adding the DMN ROIs to the atlas
        
        for i = 1:size(rois_dmntpn,1)
            comma = find(rois_dmntpn{i,1}==',');
            rois_dmntpn{i,1} = rois_dmntpn{i,1}(1:comma(1)-1);
            if i == 18 || i == 36
                rois_dmntpn{i,1} = 'cerebellum';
            end
        end
        % Fixing the labels        
        save('8-dmn_tpn_rois','rois_dmntpn');
        
    end
    % Saving all the data appropriately. It's a bit more complicated than
    % that should be since I am creating rois_consolidated, which I will
    % talk about a little bit later, and is quite central to understanding
    % the methods.

end

% ----------------------------------------------------------------------- %

% In addition to comapring the regions within the DMN and TPN for both
% groups, I also want to know if the level of anti-correlation between the
% DMN and TPN is different when compared between control individuals and
% those with ADHD. To find out, I am going to get the DMN and TPN
% timecourses from each scan, calculate the anti-correlation between them,
% concatenate the numbers for each group, and test for differences.

for dmn_tpn_anticorrelation = 1:2
    
    cd ([dir_groups{dmn_tpn_anticorrelation},'/resources'])
    load('networks.mat'), cd ..
    fprintf(['\n\n',whatsthetime,' - Acquiring DMN/TPN AC values ', ...
        'from group ',num2str(figure1c),'\n'])
    % Loading the DMN and TPN masks, which will be necessary 

    dmnsignals = [];
    tpnsignals = [];
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D.mat']);
        [x,y] = get_roi_tc(dmn,func,tpn);
        dmnsignals = cat(1,dmnsignals,x);
        tpnsignals = cat(1,tpnsignals,y);
        cd ..
    end
    % Getting the DMN and TPN timecourses from all scans
    
    correlations = zeros(length(subjects),1);
    for i = 1:length(subjects)
        x = corrcoef(dmnsignals(i,:),tpnsignals(i,:));
        correlations(i) = x(2);
    end
    % Calculating the correlation between DMN and TPN timecourses for each
    % scan/subject individually
    
    cd (dir_results)
    if dmn_tpn_anticorrelation == 1
        dmnsignals_g1 = dmnsignals;
        tpnsignals_g1 = tpnsignals;
        corrs_g1 = correlations;
        % Saving the results from group 1 for now
    else
        dmnsignals_g2 = dmnsignals;
        tpnsignals_g2 = tpnsignals;
        corrs_g2 = correlations ;
        % Saving the results from group 2 for now 
        
        stats = cell(7,2);
        stats{1,1} = 'p-value';
        stats{2,1} = 'corr-avg-g1';
        stats{3,1} = 'corr-avg-g2';
        stats{4,1} = 'corr-med-g1';
        stats{5,1} = 'corr-med-g2';
        stats{6,1} = 'corr-std-g1';
        stats{7,1} = 'corr-std-g2';
        stats{1,2} = ranksum(corrs_g1,corrs_g2);
        stats{2,2} = mean(corrs_g1);
        stats{3,2} = mean(corrs_g2);
        stats{4,2} = median(corrs_g1);
        stats{5,2} = median(corrs_g2);
        stats{6,2} = std(corrs_g1);
        stats{7,2} = std(corrs_g2);    
        % Conducting statistics and saving the results 
        
        save('7-dmn_tpn_anticorrelationn','dmnsignals_g1',  ...
            'tpnsignals_g1','dmnsignals_g2','tpnsignals_g2', ...
            'corrs_g1','corrs_g2','stats')
        % Saving all the variables from this step of methods 

        f = figure;
        a = axes;
        h = histogram(corrs_g1);
        h.LineWidth = 3;
        h.BinWidth = 0.05;
        set(f,'color','white')
        set(a,'ylim',[0 30],'xlim',[-1.1 0],'box','on','tickdir','out')
        title('Control')
        f = figure;
        a = axes;
        h = histogram(corrs_g2);
        h.LineWidth = 3;
        h.BinWidth = 0.05;
        set(f,'color','white')
        set(a,'ylim',[0 40],'xlim',[-1.1 0],'box','on','tickdir','out')
        title('ADHD')    
        % Plotting figures that will become part of Figure  1c

    end

end

% ----------------------------------------------------------------------- %

% Lastly, I want to know if the regions in the DMN mask in the control grup
% are more correlated with the PCC than the DMN mask regions in the ADHD
% group. I want to then conduct the same analysis for the TPN. My
% hypothesis is that the correlation of DMN regions with PCC will be
% stronger in control and that anticorrelation of TPN regions with PCC will
% be stronger in control than ADHD. Let's see.

for dmn_tpn_pcc_corr = 1:1
    cd (dir_results)
    load('5-dmn_tpn_info_group1.mat')
    dmnrois1 = dmn_rois;
    tpnrois1 = tpn_rois;
    load('6-dmn_tpn_info_group2.mat')
    dmnrois2 = dmn_rois;
    tpnrois2 = tpn_rois;
    x = cell2mat(dmnrois1(:,5));
    y = cell2mat(dmnrois2(:,5));
    xu = mean(x);
    xs = std(x);
    yu = mean(y);
    ys = std(y);
    [h,p] = ttest2(x,y);
    fprintf(['\n\nDMN ROIs FC w/PCC\n\nControl mean = ', num2str(xu), ...
        '\n\nControl stdv = ',num2str(xs),'\n\nADHD mean = ', ...
        num2str(yu),'\nADHD stdv = ',num2str(ys), '\n\nStatistics\n', ...
        'h = ',num2str(h),'\np = ',num2str(p),'\n\n'])
    x = cell2mat(tpnrois1(:,5));
    y = cell2mat(tpnrois2(:,5));
    xu = mean(x);
    xs = std(x);
    yu = mean(y);
    ys = std(y);
    [h,p] = ttest2(x,y);
    fprintf(['\n\nTPN ROIs FC w/PCC\n\nControl mean = ', num2str(xu), ...
        '\nControl stdv = ',num2str(xs),'\n\nADHD mean = ', ...
        num2str(yu),'\nADHD stdv = ',num2str(ys), '\n\nStatistics\n', ...
        'h = ',num2str(h),'\np = ',num2str(p),'\n\n'])    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%                 Spatial Comparison of QPPs across groups                %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Next, I am going to run the pattern-finding algorithm on the data. I will
% be concatenating 30 scans from each group like I did to acquire DMN and
% TPN masks and running the pattern-finding algorithm a hundred times and
% picking the outcome that shows the closest thing to a DMN > TPN switch
% instead of the same pattern in a different phase. I'm also going to
% create figures showing the DMN and TPN timecourses within the QPPs for
% each group and also how different they are from one another.

for acquiring_qpps = 1:2
    
    cd (groups{acquiring_qpps})
    fprintf(['\nAcquiring QPPs from group ',num2str(acquiring_qpps),'\n'])

    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects_to_use)
        cd (subjects{subjects_to_use(i)})
        if i == 1, load('t1_reorient_brain_regatlas_gm_2D.mat'), end
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D.mat']);
        f = cat(3,f,func);
        cd ..
    end
    % Concatenating the scans that the pattern finding algorithm will be
    % run on. 25 has previously been shown (Abbas et al. 2018) to be
    % sufficient for acquiring a representative QPP. Here, I am using 30.
    
    stps = [2281; 2988];
    [qpp,scv] = find_pattern_majeed2011(f,gm,10,stps(acquiring_qpps));
    cd resources
    save('qpp','qpp','scv')
    % Running the pattern findinng algorithm and saving the resulting QPP
    % in the resources folder. This isn't how the QPPs were actually
    % acquired but these starting timepoints will help you recreate the
    % results.
    
    cd (dir_results)
    [dmnsignal,tpnsignal] = get_roi_tc(dmn,qpp,tpn);
    if acquiring_qpps == 1
        qpp_brain = qpp;
        qpp_brain(qpp_brain==0) = NaN;
        thresh_upper = prctile(qpp_brain(:),90);
        thresh_lower = prctile(qpp_brain(:),10);
        qpp_upper = qpp .* (qpp > thresh_upper);
        qpp_lower = qpp .* (qpp < thresh_lower);
        qpp = qpp_upper + qpp_lower;
        unslice_save(qpp,[91 109 91 10],'9-qpp_group1.nii.gz', ...
            [2 2 2],[45 65 37]);
        % Saving the QPP as .nii.gz so that I can use it in Figure 2
        % Only keeping the top and bottom 10% of activation
        dmnsignal_g1 = dmnsignal;
        tpnsignal_g1 = tpnsignal;
        % Keeping the DMN/TPN timecourses to be saved with ADHD group's
    else
        qpp_brain = qpp;
        qpp_brain(qpp_brain==0) = NaN;
        thresh_upper = prctile(qpp_brain(:),90);
        thresh_lower = prctile(qpp_brain(:),10);
        qpp_upper = qpp .* (qpp > thresh_upper);
        qpp_lower = qpp .* (qpp < thresh_lower);
        qpp = qpp_upper + qpp_lower;
        unslice_save(qpp,[91 109 91 1],'10-qpp_group2.nii.gz', ...
            [2 2 2],[45 65 37]);
        % Saving the QPP as .nii.gz so that I can use it in Figure 2
        % Only keeping the top and bottom 10% of activation
        dmnsignal_g2 = dmnsignal;
        tpnsignal_g2 = tpnsignal;
        time = 2:2:(length(dmnsignal)*2);
        save('11-qpp_dmn_tpn_tcs','dmnsignal_g1','dmnsignal_g2',  ...
            'tpnsignal_g1','tpnsignal_g2','time');
        
        figure
        hold on
        plot(time,dmnsignal_g1,'linewidth',5);
        plot(time,tpnsignal_g1,'linewidth',5);
        set(gca,'ylim',[-0.75 0.75],'xlim',[2 20],'xtick',[],'ytick',[])
        % Plotting the DMN and TPN timecourses for the Control QPP

        figure
        hold on
        plot(time,dmnsignal_g2,'linewidth',5);
        plot(time,tpnsignal_g2,'linewidth',5);
        set(gca,'ylim',[-0.75 0.75],'xlim',[2 20],'xtick',[],'ytick',[])
        % Plotting the DMN and TPN timecourses for the ADHD QPP

        figure
        hold on
        dmntpn_diff_g1 = zeros(length(dmnsignal_g1),1);
        dmntpn_diff_g2 = zeros(length(dmnsignal_g2),1);
        for i = 1:length(dmnsignal_g1)
            dmntpn_diff_g1(i) = (dmnsignal_g1(i)-tpnsignal_g1(i))^2;
            dmntpn_diff_g2(i) = (dmnsignal_g2(i)-tpnsignal_g2(i))^2;
        end
        plot(time,dmntpn_diff_g1,'linewidth',5);
        plot(time,dmntpn_diff_g2,'linewidth',5);
        set(gca,'ylim',[0 1.5],'xlim',[2 20],'xtick',[],'ytick',[])
        % Plotting the square of the difference between the DMN timecourse 
        % and the TPN timecourse in the Control and ADHD QPPs
        
    end
    
end

% ----------------------------------------------------------------------- %

% Now that I have representative QPPs from  each group, I want to
% compare them to see if there are any differences. I am again going to
% use the Brainnetome atlas here. Instead of comparing the spatial
% pattern voxel by voxel, I am going to do for the 273 ROIs in the atlas.
% The results of this analysis will also end up in Supplementary Table 2.

for comparing_qpps = 1:1
    
    cd (dir_group1), cd resources
    load('qpp.mat','qpp')
    qpp_g1 = qpp; clear qpp
    cd (dir_group2), cd resources
    load('qpp.mat','qpp')
    qpp_g2 = qpp; clear qpp
    cd (dir_atlas)
    load('roi_atlas.mat')
    load('labels.mat')
    % Loading all the files I'm going to need 
    
    qpp_group1_tcs = zeros(length(rois),size(qpp_g1,3));
    qpp_group2_tcs = zeros(length(rois),size(qpp_g2,3));
    tc_corrs = zeros(length(rois),1);
    for i = 1:length(rois)
        qpp_group1_tcs(i,:) = get_roi_tc(rois{i},qpp_g1);
        qpp_group2_tcs(i,:) = get_roi_tc(rois{i},qpp_g2);
        x = corrcoef(qpp_group1_tcs(i,:),qpp_group2_tcs(i,:));
        tc_corrs(i) = x(2);
    end
    % Calculating correlation for all ROIs 

    corrmap = zeros(size(qpp_g1,1),size(qpp_g1,2));
    for i = 1:length(rois)
        if tc_corrs(i) < -0.3 || tc_corrs(i) > 0.3
            corrmap(rois{i}==1) = tc_corrs(i);
        end
    end
    cd (dir_results)
    unslice_save(corrmap,[91 109 91],'12-qpp_comparison.nii.gz',[2 2 2],...
        [45 65 37]);
    % Creating and saving the correlation map that will end up in the
    % figure and will nicely show the similarities and differences in one
    % image
    
    load('5-dmn_tpn_info_group1.mat');
    dmnrois1 = dmn_rois;
    tpnrois1 = tpn_rois;
    load('6-dmn_tpn_info_group2.mat');
    dmnrois2 = dmn_rois;
    tpnrois2 = tpn_rois;
    roi_corrs = cell(273,6);
    for i = 1:273
        roi_corrs{i,1} = labels{i};
        for j = 1:length(dmnrois1)
            if strcmp(roi_corrs{i,1},dmnrois1{j,2}) == 1
                roi_corrs{i,2} = 1;
            end
        end
        for j = 1:length(tpnrois1)
            if strcmp(roi_corrs{i,1},tpnrois1{j,2}) == 1
                roi_corrs{i,4} = 1;
            end
        end
        for j = 1:length(dmnrois2)
            if strcmp(roi_corrs{i,1},dmnrois2{j,2}) == 1
                roi_corrs{i,3} = 1;
            end
        end
        for j = 1:length(tpnrois2)
            if strcmp(roi_corrs{i,1},tpnrois2{j,2}) == 1
                roi_corrs{i,5} = 1;
            end
        end
        roi_corrs{i,6} = tc_corrs(i);
    end
    % This is the table that will make up Supplementary Table 3
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%                Temporal Comparison of QPPs across groups                %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% To test the effect that the QPPs have on functional connectivity, I am
% going to regress the QPPs from the functional scans and calculate
% functional connectivity strength before and after I do the regression. I
% already know that functional connectivity in the DMN and TPN is partially
% dependent on the presence of the QPPs. However, it will be interesting to
% see how the regression affects the control and ADHD groups differently.

for regressing_qpps = 1:2
    fprintf(['\nRegressing QPPs from group ',num2str(regressing_qpps), ...
        '\n\n'])
    cd (dir_group1), cd resources, load('qpp.mat','qpp'), qpp_g1 = qpp;
    cd (dir_group2), cd resources, load('qpp.mat','qpp'), qpp_g2 = qpp;
    % Loading the QPPs from both groups

    cd (dir_groups{regressing_qpps})
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        fprintf([whatsthetime,' - ',subjects{i},'\n'])
        cd (subjects{i})
        is_done = dir('*qppreg.mat');
        if isempty(is_done) == 1
            load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
                '_3D.mat']);
            scv_g1 = get_swc(func,qpp_g1);
            scv_g2 = get_swc(func,qpp_g2);
            func_g1r = regress_pattern(func,qpp_g1,scv_g1);
            func_g2r = regress_pattern(func,qpp_g2,scv_g2);
            scv_g1r = get_swc(func_g1r,qpp_g1);
            scv_g2r = get_swc(func_g2r,qpp_g2);
            save(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
                '_3D_qppreg'],'func_g1r','func_g2r','scv_g1','scv_g1r', ...
                'scv_g2','scv_g2r');
            cd ..
        else
            cd ..
            continue
        end
    end
    % Regressing both QPPs from each scan in each group
end

% ----------------------------------------------------------------------- %

% To compare the temporal pattern of the QPPs, I need to look at the
% sliding correlation vectors. The distribution of values in those vectors
% tell me how often the QPP is occurring and at what correlation strength.
% So, I'm now going to load all of the sliding correlation vectors and save
% them so they can be compared cumulatively.

for gathering_scvs = 1:2
    
    if gathering_scvs == 1
        control_qpp_in_control = [];
        adhd_qpp_in_control = [];
        control_qpp_in_control_reg = [];
        adhd_qpp_in_control_reg = [];
        control_qpp_in_adhd = [];
        adhd_qpp_in_adhd = [];
        control_qpp_in_adhd_reg = [];
        adhd_qpp_in_adhd_reg = [];
    end
    % Predefining all the vectors that will hold the cumulative sliding
    % correlation vectors for both the control and ADHD QPPs
    
    cd (dir_groups{gathering_scvs})
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D_qppreg'],'scv_g1','scv_g2','scv_g1r','scv_g2r')
        if group == 1
            control_qpp_in_control = cat(1,control_qpp_in_control,scv_g1);
            adhd_qpp_in_control = cat(1,adhd_qpp_in_control,scv_g2);
            control_qpp_in_control_reg = cat(1,control_qpp_in_control_reg,scv_g1r);
            adhd_qpp_in_control_reg = cat(1,adhd_qpp_in_control_reg,scv_g2r);
        else
            control_qpp_in_adhd = cat(1,control_qpp_in_adhd,scv_g1);
            adhd_qpp_in_adhd = cat(1,adhd_qpp_in_adhd,scv_g2);
            control_qpp_in_adhd_reg = cat(1,control_qpp_in_adhd_reg,scv_g1r);
            adhd_qpp_in_adhd_reg = cat(1,adhd_qpp_in_adhd_reg,scv_g2r);
        end
        cd ..
    end
    % Concatenating the sliding correlation vectors from each scan to
    % fill in the variables I predefined before this step
    
    if gathering_scvs == 2
        cd (dir_results)
        save('13-scv_vectors','control_qpp_in_control', ...
            'adhd_qpp_in_control','control_qpp_in_control_reg', ...
            'adhd_qpp_in_control_reg','control_qpp_in_adhd', ...
            'adhd_qpp_in_adhd','control_qpp_in_adhd_reg', ...
            'adhd_qpp_in_adhd_reg');
    end
    % Saving the correlation vectors 
    
end

% ----------------------------------------------------------------------- %

% Figure 3a will show examples of what the sliding correlation vectors look
% like - both before and after QPP regression. Figure 3b will count the
% mean height of peaks greater than 0.1 in correlation for control and ADHD
% before regression and after. Figure 3c will count the mean frequency of
% peaks greater than 0.1 in correlation for control and ADHD befor
% regression and after. Figure 3c will count the mean frequency of peaks
% greater than 0.1 in correlation for control and ADHD befor regression and
% after.

for figure_3  = 1:1
    
    cd (dir_results)
    load('13-scv_vectors.mat')
    
    % ------------------------------ (a) -------------------------------- % 

    figure 
    hold on
    plot(control_qpp_in_control,'linewidth',3)
    plot(control_qpp_in_control_reg,'linewidth',3)
    set(gca,'box','on','tickdir','out','ylim',[-0.75 0.75],'xlim',[1 600])

    figure 
    hold on
    plot(adhd_qpp_in_control,'linewidth',3)
    plot(adhd_qpp_in_control_reg,'linewidth',3)
    set(gca,'box','on','tickdir','out','ylim',[-0.75 0.75],'xlim',[1 600])
    
    % ------------------------------ (b) -------------------------------- % 


    figure 
    hold on
    plot(control_qpp_in_adhd,'linewidth',3)
    plot(control_qpp_in_adhd_reg,'linewidth',3)
    set(gca,'box','on','tickdir','out','ylim',[-0.75 0.75],'xlim',[1 600])

    figure 
    hold on
    plot(adhd_qpp_in_adhd,'linewidth',3)
    plot(adhd_qpp_in_adhd_reg,'linewidth',3)
    set(gca,'box','on','tickdir','out','ylim',[-0.75 0.75],'xlim',[1 600])
    
    % ------------------------------ (c) -------------------------------- % 
    
    pk_height_means_g1 = zeros(1,2);
    % column 1 - strength of control QPP in control scans 
    % column 2 - strength of control QPP in adhd scans
    pk_height_means_g1r = zeros(1,2);
    % column 1 - strength of control QPP in QPP regressed control scans 
    % column 2 - strength of control QPP in QPP regressed adhd scans 
    pk_height_means_g2 = zeros(1,2);
    % column 1 - strength of adhd QPP in control scans
    % column 2 - strength of adhd QPP in adhd scans
    pk_height_means_g2r = zeros(1,2);
    % column 1 - strength of adhd QPP in QPP regressed control scans 
    % column 2 - strength of adhd QPP in QPP regressed adhd scans
    
    cd (dir_group1)
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        fprintf([subjects{i},'\n'])
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D_qppreg'],'scv_g1','scv_g2','scv_g1r','scv_g2r')
        [pks,~] = findpeaks(scv_g1,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g1(i,1) = mean(pks,'omitnan');
        [pks,~] = findpeaks(scv_g1r,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g1r(i,1) = mean(pks,'omitnan');
        [pks,~] = findpeaks(scv_g2,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g2(i,1) = mean(pks,'omitnan');
        [pks,~] = findpeaks(scv_g2r,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g2r(i,1) = mean(pks,'omitnan');
        cd ..
    end
    % Getting the mean peak height for all peaks greater than 0.1 for each
    % sliding correlation vector
    
    cd (dir_group2)
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        fprintf([subjects{i},'\n'])
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D_qppreg'],'scv_g1','scv_g2','scv_g1r','scv_g2r')
        [pks,~] = findpeaks(scv_g1,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g1(i,2) = mean(pks,'omitnan');
        [pks,~] = findpeaks(scv_g1r,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g1r(i,2) = mean(pks,'omitnan');
        [pks,~] = findpeaks(scv_g2,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g2(i,2) = mean(pks,'omitnan');
        [pks,~] = findpeaks(scv_g2r,'minpeakheight',0.1);
        if isempty(pks) == 1, pks = NaN; end
        pk_height_means_g2r(i,2) = mean(pks,'omitnan');
        cd ..
    end
    % Getting the mean peak height for all peaks greater than 0.1 for each
    % sliding correlation vector
    
    [h1,p1] = ttest2(pk_height_means_g1(:,1),pk_height_means_g1r(:,1));
    [h2,p2] = ttest2(pk_height_means_g2(:,1),pk_height_means_g2r(:,1));
    [h3,p3] = ttest2(pk_height_means_g1(:,1),pk_height_means_g2(:,1));
    [h4,p4] = ttest2(pk_height_means_g1(:,2),pk_height_means_g1r(:,2));
    [h5,p5] = ttest2(pk_height_means_g2(:,2),pk_height_means_g2r(:,2));
    [h6,p6] = ttest2(pk_height_means_g1(:,2),pk_height_means_g2(:,2));
    [h7,p7] = ttest2(pk_height_means_g1(:,1),pk_height_means_g2(:,2));
    [h8,p8] = ttest2(pk_height_means_g1(:,2),pk_height_means_g2(:,1));
    cd (dir_results)
    save('14-qpp_strength','pk_height_means_g1', ...
        'pk_height_means_g1r','pk_height_means_g2', ...
        'pk_height_means_g2r','h1','h2','h3','h4','h5','h6','h7', ...
        'p1','p2','p3','p4','p5','p6','p7','p8','h8');
    % Conducting statistics and saving the results
    
    % ------------------------------ (d) -------------------------------- % 
    
    pk_freq_means_g1 = zeros(1,2);
    % column 1 - freq of cont QPP in cont scans 
    % column 2 - freq of cont QPP in adhd scans
    pk_freq_means_g1r = zeros(1,2);
    % column 1 - freq of cont QPP in cont QPP regressed cont scans 
    % column 2 - freq of cont QPP in cont QPP regressed adhd scans 
    pk_freq_means_g2 = zeros(1,2);
    % column 1 - freq of adhd QPP in cont scans`
    % column 2 - freq of adhd QPP in adhd scans
    pk_freq_means_g2r = zeros(1,2);
    % column 1 - freq of adhd QPP in adhd QPP regressed cont scans 
    % column 2 - freq of adhd QPP in adhd QPP regressed adhd scans
    
    cd (dir_group1)
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D_qppreg'],'scv_g1','scv_g2','scv_g1r','scv_g2r')        
        [pks,~] = findpeaks(scv_g1,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g1(i,1) = NaN;
        else
           pk_freq_means_g1(i,1) = mean((length(pks)/5.33),'omitnan');
        end
        [pks,~] = findpeaks(scv_g1r,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g1r(i,1) = NaN;
        else
          pk_freq_means_g1r(i,1) = mean((length(pks)/5.33),'omitnan');
        end
        [pks,~] = findpeaks(scv_g2,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g2(i,1) = NaN;
        else
            pk_freq_means_g2(i,1) = mean((length(pks)/5.33),'omitnan');
        end
        [pks,locs] = findpeaks(scv_g2r,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g2r(i,1) = NaN;
        else
           pk_freq_means_g2r(i,1) = mean((length(pks)/5.33),'omitnan');
        end
        cd ..
    end
    % Calculating the mean peak distances
    
    cd (dir_group2)
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D_qppreg'],'scv_g1','scv_g2','scv_g1r','scv_g2r')        
        [pks,~] = findpeaks(scv_g1,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g1(i,2) = NaN;
        else
           pk_freq_means_g1(i,2) = mean((length(pks)/5.33),'omitnan');
        end
        [pks,~] = findpeaks(scv_g1r,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g1r(i,2) = NaN;
        else
          pk_freq_means_g1r(i,2) = mean((length(pks)/5.33),'omitnan');
        end
        [pks,~] = findpeaks(scv_g2,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g2(i,2) = NaN;
        else
            pk_freq_means_g2(i,2) = mean((length(pks)/5.33),'omitnan');
        end
        [pks,locs] = findpeaks(scv_g2r,'minpeakheight',0.1);
        if isempty(pks) == 1
            pk_freq_means_g2r(i,2) = NaN;
        else
           pk_freq_means_g2r(i,2) = mean((length(pks)/5.33),'omitnan');
        end
        cd ..
    end
    % Calculating the mean peak distances  
    
    [h1,p1] = ttest2(pk_freq_means_g1(:,1),pk_freq_means_g1r(:,1));
    [h2,p2] = ttest2(pk_freq_means_g2(:,1),pk_freq_means_g2r(:,1));
    [h3,p3] = ttest2(pk_freq_means_g1(:,1),pk_freq_means_g2(:,1));
    [h4,p4] = ttest2(pk_freq_means_g1(:,2),pk_freq_means_g1r(:,2));
    [h5,p5] = ttest2(pk_freq_means_g2(:,2),pk_freq_means_g2r(:,2));
    [h6,p6] = ttest2(pk_freq_means_g1(:,2),pk_freq_means_g2(:,2));
    [h7,p7] = ttest2(pk_freq_means_g1(:,1),pk_freq_means_g2(:,2));
    [h8,p8] = ttest2(pk_freq_means_g1(:,2),pk_freq_means_g2(:,1));
    cd (dir_results)
    save('15-qpp_frequency','pk_freq_means_g1', ...
        'pk_freq_means_g1r','pk_freq_means_g2', ...
        'pk_freq_means_g2r','h1','h2','h3','h4','h5','h6','h7', ...
        'h8','p1','p2','p3','p4','p5','p6','p7','p8');
    % Saving the results
    
    % ------------------------------ (e) -------------------------------- % 
    
    f1 = figure;
    f1.Color = 'white';
    ax = axes;
    hold on
    h1 = histogram(control_qpp_in_control);
    h1.LineWidth = 3;
    h1.BinLimits = [-0.75 0.75];
    h1.BinWidth = 0.0375;
    h1.FaceColor = [0.3 0.3 1];
    h2 = histogram(control_qpp_in_control_reg);
    h2.LineWidth = 3;
    h2.BinLimits = [-0.75 0.75];
    h2.BinWidth = 0.0375;
    h2.FaceColor = [1 0.3 0.3];
    ax.XLim = [-0.75 0.75];
    ax.YLim = [1 5000];
    ax.Box = 'on';
    ax.LineWidth = 3;
    ax.TickDir = 'out';
    xlabel('Sliding Correlation');
    ylabel('Frequency (Hz)');
    
    f2 = figure;
    f2.Color = 'white';
    ax = axes;
    hold on
    h1 = histogram(control_qpp_in_adhd);
    h1.LineWidth = 3;
    h1.BinLimits = [-0.75 0.75];
    h1.BinWidth = 0.0375;
    h1.FaceColor = [0.3 0.3 1];
    h2 = histogram(control_qpp_in_adhd_reg);
    h2.LineWidth = 3;
    h2.BinLimits = [-0.75 0.75];
    h2.BinWidth = 0.0375;
    h2.FaceColor = [1 0.3 0.3];
    ax.XLim = [-0.75 0.75];
    ax.YLim = [1 5000];
    ax.Box = 'on';
    ax.LineWidth = 3;
    ax.TickDir = 'out';
    xlabel('Sliding Correlation');
    ylabel('Frequency (Hz)');
    
    % ------------------------------ (f) -------------------------------- % 
    
    f3 = figure;
    f3.Color = 'white';
    ax = axes;
    hold on
    h1 = histogram(adhd_qpp_in_control);
    h1.LineWidth = 3;
    h1.BinLimits = [-0.75 0.75];
    h1.BinWidth = 0.0375;
    h1.FaceColor = [0.3 0.3 1];
    h2 = histogram(adhd_qpp_in_control_reg);
    h2.LineWidth = 3;
    h2.BinLimits = [-0.75 0.75];
    h2.BinWidth = 0.0375;
    h2.FaceColor = [1 0.3 0.3];
    ax.XLim = [-0.75 0.75];
    ax.YLim = [1 5000];
    ax.Box = 'on';
    ax.LineWidth = 3;
    ax.TickDir = 'out';
    xlabel('Sliding Correlation');
    ylabel('Frequency (Hz)');
    
    f4 = figure;
    f4.Color = 'white';
    ax = axes;
    hold on
    h1 = histogram(adhd_qpp_in_adhd);
    h1.LineWidth = 3;
    h1.BinLimits = [-0.75 0.75];
    h1.BinWidth = 0.0375;
    h1.FaceColor = [0.3 0.3 1];
    h2 = histogram(adhd_qpp_in_adhd_reg);
    h2.LineWidth = 3;
    h2.BinLimits = [-0.75 0.75];
    h2.BinWidth = 0.0375;
    h2.FaceColor = [1 0.3 0.3];
    ax.XLim = [-0.75 0.75];
    ax.YLim = [1 5000];
    ax.Box = 'on';
    ax.LineWidth = 3;
    ax.TickDir = 'out';
    xlabel('Sliding Correlation');
    ylabel('Frequency (Hz)');
    
    [h1,p1] = kstest2(control_qpp_in_control,control_qpp_in_control_reg);
    [h2,p2] = kstest2(control_qpp_in_adhd,control_qpp_in_adhd_reg);
    [h3,p3] = kstest2(control_qpp_in_control,control_qpp_in_adhd);
    [h4,p4] = kstest2(adhd_qpp_in_control,adhd_qpp_in_control_reg);
    [h5,p5] = kstest2(adhd_qpp_in_adhd,adhd_qpp_in_adhd_reg);
    [h6,p6] = kstest2(adhd_qpp_in_control,adhd_qpp_in_adhd);
    [h7,p7] = kstest2(control_qpp_in_control,adhd_qpp_in_adhd);
    [h8,p8] = kstest2(control_qpp_in_adhd,adhd_qpp_in_control);
    save('16-hist_stats','h1','h2','h3','h4','h5','h6','h7','h8', ...
        'p1','p2','p3','p4','p5','p6','p7','p8')

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%                  Comparison of Funnctional Connectivity                 %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% When acquiring DMN and TPN masks, I saved an updated ROI atlas that is
% different from the Brainnetome atlas and specific to this study. I kept
% only the ROIs that had voxels in them that belonged either to the DMN or
% TPN. Within each of those ROIs, I discarded the voxels in the mask that
% were not in the DMN or TPN. Finally, I consolidated the ROIs by the
% region hierarchy laid out in the Brainnetome atlas. Now, I'm going to
% calculate the functional connectivity matrices using this new atlas.

for calculating_fc_matrices = 1:2
    
    fprintf(['\n\nCalculating FC matrices for group ',  ...
        num2str(calculating_fc_matrices),'\n'])
    cd (dir_results), load('8-dmn_tpn_rois.mat');
    rois = rois_consolidated;
    num_rois = length(rois);
    % Loading the DMN TPN ROI atlas
    
    if calculating_fc_matrices == 1
        fc_g1 = [];
        fc_g1_g1r = [];
        fc_g1_g2r = [];
        fc_g2 = [];
        fc_g2_g1r = [];
        fc_g2_g2r = [];
    end    
    %  Predefining the matrices that will hold all the functional
    %  connectivity matrices
    
    count = 0;    
    cd (dir_groups{calculating_fc_matrices})
    subjects = dir('subject*');
    subjects = {subjects.name}';
    for i = 1:length(subjects)
        cd (subjects{i})
        fprintf([whatsthetime,' - ',subjects{i},'\n']);
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D.mat']);
        load(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D_qppreg.mat'],'func_g1r','func_g2r');
        roi_tcs = zeros(num_rois,size(func,3));
        roi_tcs_g1r = zeros(num_rois,size(func,3));
        roi_tcs_g2r = zeros(num_rois,size(func,3));
        parfor k = 1:num_rois
            
            roi_tcs(k,:) = get_roi_tc(rois{k},func);
            roi_tcs_g1r(k,:) = get_roi_tc(rois{k},func_g1r);
            roi_tcs_g2r(k,:) = get_roi_tc(rois{k},func_g2r);
            
        end
        fc_matrix = ones(num_rois,num_rois);
        fc_matrix_g1r = ones(num_rois,num_rois);
        fc_matrix_g2r = ones(num_rois,num_rois);
        for k = 1:num_rois
            
            if k ~= num_rois
                timesaver = k + 1;
            else
                continue
            end
            
            tc = roi_tcs(k,:);
            tc_g1r = roi_tcs_g1r(k,:);
            tc_g2r = roi_tcs_g2r(k,:);
            
            parfor l = timesaver:num_rois
                
                fc_strength = corrcoef(tc,roi_tcs(l,:));
                fc_matrix(k,l) = fc_strength(2);
                
                fc_strength = corrcoef(tc_g1r,roi_tcs_g1r(l,:));
                fc_matrix_g1r(k,l) = fc_strength(2);
                
                fc_strength = corrcoef(tc_g2r,roi_tcs_g2r(l,:));
                fc_matrix_g2r(k,l) = fc_strength(2);
                
            end
            
            fc_matrix(:,k) = fc_matrix(k,:)';
            fc_matrix_g1r(:,k) = fc_matrix_g1r(k,:)';
            fc_matrix_g2r(:,k) = fc_matrix_g2r(k,:)';
            
        end
        if calculating_fc_matrices == 1
            fc_g1 = cat(3,fc_g1,fc_matrix);
            fc_g1_g1r = cat(3,fc_g1_g1r,fc_matrix_g1r);
            fc_g1_g2r = cat(3,fc_g1_g2r,fc_matrix_g2r);
        else
            fc_g2 = cat(3,fc_g2,fc_matrix);
            fc_g2_g1r = cat(3,fc_g2_g1r,fc_matrix_g1r);
            fc_g2_g2r = cat(3,fc_g2_g2r,fc_matrix_g2r);
        end
        save(['f01_reorient_trim_stc_mc_regatlas_sm_fil_gsr_wmcsfr_zsc',...
            '_3D_qppreg_fc_dmntpn.mat'],'fc_matrix','fc_matrix_g1r', ...
            'fc_matrix_g2r');
        cd ..
    end
    

    if calculating_fc_matrices == 2
        cd (dir_results)
        for i = 1:size(fc_g1,1)
            for j = 1:size(fc_g1,2)
                for k = 1:size(fc_g1,3)
                    fc_g1(i,j,k) = 0.5 * ...
                        (log(1+fc_g1(i,j,k))- ...
                        log(1-fc_g1(i,j,k))); 
                    fc_g1_g1r(i,j,k) = 0.5 * ...
                        (log(1+fc_g1_g1r(i,j,k))- ...
                        log(1-fc_g1_g1r(i,j,k)));
                    fc_g1_g2r(i,j,k) = 0.5 * ...
                        (log(1+fc_g1_g2r(i,j,k))- ...
                        log(1-fc_g1_g2r(i,j,k)));
                    fc_g2(i,j,k) = 0.5 * ...
                        (log(1+fc_g2(i,j,k))- ...
                        log(1-fc_g2(i,j,k))); 
                    fc_g2_g1r(i,j,k) = 0.5 * ...
                        (log(1+fc_g2_g1r(i,j,k))- ...
                        log(1-fc_g2_g1r(i,j,k)));
                    fc_g2_g2r(i,j,k) = 0.5 * ...
                        (log(1+fc_g2_g2r(i,j,k))- ...
                        log(1-fc_g2_g2r(i,j,k))); %#ok<*SAGROW>
                end
            end
        end
        cd (dir_results)
        save('17-fc_matrices','fc_g1','fc_g1_g1r', ...
            'fc_g1_g2r','fc_g2','fc_g2_g1r', ...
            'fc_g2_g2r');
    end
    % Doing a fischer z-transform and saving the matrices    
end

% ----------------------------------------------------------------------- %

% Figure 4a is going to show the mean functional connectivity in the
% control scans and ADHD scans. Figure 4b is going to show the significant
% differences between the two groups before their native QPPs were
% regressed and after their native QPPs were regressed. Figure 4c is
% going to show the effect of the regression of the control QPP and ADHD
% QPP on the control scans. Figure 4d is going to show the effect of the
% regression of the control QPP and ADHD QPP on the ADHD scans.

for figure4 = 1:1
    
    cd (dir_results)
    load('8-dmn_tpn_rois.mat');
    load('17-fc_matrices.mat');
    num_rois = size(fc_g1,1);
    

    % ------------------------------ (a) -------------------------------- % 
    
    % This is the part of the figure that will be showing the mean
    % functional connectivity in the control group and ADHD group. No
    % statistics yet; just displaying the mean of all 106 functional scans
    % from each group. The bottom left shows Control and top right shows
    % ADHD.
    
    fc_g1_mean = mean(fc_g1,3);
    fc_g2_mean = mean(fc_g2,3);
    fcg1mean_and_fcg2mean = combine_fc_matrices(fc_g1_mean,fc_g2_mean);
    
    figure
    imagesc(fcg1mean_and_fcg2mean)
    colormap jet
    set(gca,'xtick',1:36,'ytick',1:36,'tickdir','out')
    set(gcf,'color','white')
    caxis([-1 1])
    xticklabels(rois_dmntpn(:,1));
    xtickangle(90);
    yticklabels(rois_dmntpn(:,1));
    
    % ------------------------------ (b) -------------------------------- % 
    
    % The next figure is a direct comparison of functional connectivity
    % between the Control and ADHD. The bottom left is a comparison of
    % functional connectivity before anything was done to the scans
    % (original scans). The top right comparison is the functional
    % connectivity differences after the native QPPs were regressed.
    
    fcg1_vs_fcg2 = compare_fc_matrices(fc_g1,fc_g2,0.01);
    fcg1g1r_vs_fcg2g2r = compare_fc_matrices(fc_g1_g1r,fc_g2_g2r,0.01);
    fcg1vsfcg2_and_fcg1g1rvsfcg2g2r = combine_fc_matrices(fcg1_vs_fcg2, ...
        fcg1g1r_vs_fcg2g2r);
    
    figure
    imagesc(fcg1vsfcg2_and_fcg1g1rvsfcg2g2r)
    colormap jet
    set(gca,'xtick',1:36,'ytick',1:36,'tickdir','out')
    set(gcf,'color','white')
    caxis([-0.3 0.3])    
    xticklabels(rois_dmntpn(:,1));
    xtickangle(90);
    yticklabels(rois_dmntpn(:,1));
    
    % ------------------------------ (c) -------------------------------- %
    
    % For the actual figure, I am plotting just what happens to functional
    % connectivity in the control functional scans after the control QPP is
    % regressed from them. In the supplementary figure, it is that in
    % addition to what happens after regression of the ADHD QPP. The latter
    % is not very relevant to the paper so i'm putting that in the
    % supplementary figure on Shella's advice.
    
    fcg1_vs_fcg1g1r = compare_fc_matrices(fc_g1,fc_g1_g1r,0.01);
    
    figure
    imagesc(fcg1_vs_fcg1g1r)
    colormap jet
    set(gca,'xtick',1:36,'ytick',1:36,'tickdir','out')
    set(gcf,'color','white')
    caxis([-0.5 0.5])    
    xticklabels(rois_dmntpn(:,1));
    xtickangle(90);
    yticklabels(rois_dmntpn(:,1));
    
    % ------- %
    
    fcg1_vs_fcg1g1r = compare_fc_matrices(fc_g1,fc_g1_g1r,0.01);
    fcg1_vs_fcg1g2r = compare_fc_matrices(fc_g1,fc_g1_g2r,0.01);
    fcg1vsfcg1g1r_and_fcg1vsfcg1g2r = combine_fc_matrices  ...
        (fcg1_vs_fcg1g1r, fcg1_vs_fcg1g2r);
    
    figure
    imagesc(fcg1vsfcg1g1r_and_fcg1vsfcg1g2r)
    colormap jet
    set(gca,'xtick',1:36,'ytick',1:36,'tickdir','out')
    set(gcf,'color','white')
    caxis([-0.5 0.5])    
    xticklabels(rois_dmntpn(:,1));
    xtickangle(90);
    yticklabels(rois_dmntpn(:,1));
    
    % ------------------------------ (d) -------------------------------- %
    
    % For the actual figure, I am plotting just what happens to functional
    % connectivity in the ADHD functional scans after the ADHD QPP is
    % regressed from them. In the supplementary figure, it is that in
    % addition to what happens after regression of the Control QPP. The 
    % latter is not very relevant to the paper so i'm putting that in the
    % supplementary figure on Shella's advice.
    
    fcg2_vs_fcg2g2r = compare_fc_matrices(fc_g2,fc_g2_g2r,0.01);
    
    figure
    imagesc(fcg2_vs_fcg2g2r)
    colormap jet
    set(gca,'xtick',1:36,'ytick',1:36,'tickdir','out')
    set(gcf,'color','white')
    caxis([-0.5 0.5])        
    xticklabels(rois_dmntpn(:,1));
    xtickangle(90);
    yticklabels(rois_dmntpn(:,1));
    
    % ------- %
    
    fcg2_vs_fcg2g1r = compare_fc_matrices(fc_g2,fc_g2_g1r,0.01);
    fcg2_vs_fcg2g2r = compare_fc_matrices(fc_g2,fc_g2_g2r,0.01);
    fcg2vsfcg2g1r_and_fcg2vsfcg2g2r = combine_fc_matrices  ...
        (fcg2_vs_fcg2g1r, fcg2_vs_fcg2g2r);
    
    figure
    imagesc(fcg2vsfcg2g1r_and_fcg2vsfcg2g2r)
    colormap jet
    set(gca,'xtick',1:36,'ytick',1:36,'tickdir','out')
    set(gcf,'color','white')
    caxis([-0.5 0.5])        
    xticklabels(rois_dmntpn(:,1));
    xtickangle(90);
    yticklabels(rois_dmntpn(:,1));
    
    % ---------------------------- extra -------------------------------- %
    
    % Here I am comparing the overall functional connectivity difference
    % between local connectivity in DMN between groups, local connectivity
    % in TPN between groups, and inter-network connectivity (or
    % anti-correlation between groups. Then I am plotting it out as
    % histograms to really show the difference.
    
    dmn_fc_g1 = fc_g1(1:18,1:18,:);
    dmn_fc_g1(dmn_fc_g1==Inf) = [];
    dmn_fc_g2 = fc_g2(1:18,1:18,:);
    dmn_fc_g2(dmn_fc_g2==Inf) = [];
    figure, hold on
    histogram(dmn_fc_g1(:)); set(gca,'xlim',[-1 1]);
    histogram(dmn_fc_g2(:)); set(gca,'xlim',[-1 1]);
    [~,p] = ttest2(dmn_fc_g1(:),dmn_fc_g2(:),0.01);
    fprintf(['\n',num2str(p),'\n'])

    tpn_fc_g1 = fc_g1(19:36,19:36,:);
    tpn_fc_g1(tpn_fc_g1==Inf) = [];
    tpn_fc_g2 = fc_g2(19:36,19:36,:);
    tpn_fc_g2(tpn_fc_g2==Inf) = [];
    figure, hold on
    histogram(tpn_fc_g1(:)), set(gca,'xlim',[-1 1])
    histogram(tpn_fc_g2(:)), set(gca,'xlim',[-1 1])
    [~,p] = ttest2(tpn_fc_g1(:),tpn_fc_g2(:),0.01);
    fprintf(['\n',num2str(p),'\n'])

    dmntpn_fc_g1 = fc_g1(1:18,19:36,:);
    dmntpn_fc_g1(dmn_fc_g1==Inf) = [];
    dmntpn_fc_g2 = fc_g2(1:18,19:36,:);
    dmntpn_fc_g2(dmn_fc_g2==Inf) = [];
    figure, hold on
    histogram(dmntpn_fc_g1(:)), set(gca,'xlim',[-1 1])
    histogram(dmntpn_fc_g2(:)), set(gca,'xlim',[-1 1])
    [~,p] = ttest2(dmntpn_fc_g1(:),dmntpn_fc_g2(:),0.01);
    fprintf(['\n',num2str(p),'\n'])    
    
    figure
    hold on
    h1 = histogram(dmn_fc_g1(:));
    h2 = histogram(dmn_fc_g2(:)); 
    
end