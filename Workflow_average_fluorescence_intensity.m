% Prerequirements (find details in "FIJI workflow for fluorescence intensity measurement" on Overleaf)
% 	- Correct background 
% 	- Correct shading (optional): no need to do this if the illumination is even 
% 	- Measure flurecence intensity of ROIs in BG and Shading corrected image 
% 
% Notes: Name all saved csv files with keywords representing the ROI types, such as 'BG', 'subnuclei', or'soma'


% Contents
%% ====================
% 1. Default settings
% clearvars -except recdata_organized alignedData_allTrials alignedData_event_list seriesData_sync grouped_event adata grouped_event_info_filtered

% Folders
PC_name = getenv('COMPUTERNAME'); 
% set folders for different situation
DataFolder = 'G:\Workspace\Inscopix_Seagate';

% FIM: fluorescence intensity measurement
if strcmp(PC_name, 'GD-AW-OFFICE')
	FIM_folder = 'D:\guoda\Documents\Workspace\Analysis\Confocal'; % office desktop
elseif strcmp(PC_name, 'LAPTOP-84IERS3H')
	FIM_folder = 'C:\Users\guoda\Documents\Workspace\Analysis\Confocal'; % laptop
end
DefaultDir.csv = FIM_folder;
DefaultDir.mat = FIM_folder;
% [FolderPathVA] = set_folder_path_ventral_approach(DataFolder,FIM_folder);


% Keywords and labels. 
% Note: Name the CSV files using the labels below to group data automatically
DilutionGroup = {'d1','d2','d3'}; 
DayGroup = {'7','14','21'};
RegionGroup = {'subnuclei','soma'};

Di_num = numel(DilutionGroup);
D_num = numel(DayGroup);
R_num = numel(RegionGroup);

group_num = Di_num*D_num*R_num;
keywords_groups = cell(1,group_num); % each cell contains a single string, such as 'D1-14-soma'
keywords_cells = cell(1,group_num); % each cell contains a cell of keywords, such as {'D1','14','soma'}
combined_csv_data = cell(1,group_num);
kgc = 1; % current keyword group counting

group_num_MonoRegion = Di_num*D_num;
keywords_groups_MonoRegion = cell(1,group_num_MonoRegion);
kgc_mr = 1;
for din = 1:Di_num
	DiStr = DilutionGroup{din};
	for dn = 1:D_num
		DStr = DayGroup{dn};
		keywords_groups_MonoRegion{kgc_mr} = sprintf('%s-%s', DiStr,DStr);
		kgc_mr = kgc_mr+1;
		for rn = 1:R_num
			RStr = RegionGroup{rn};
			keywords_groups{kgc} = sprintf('%s-%s-%s',DiStr,DStr,RStr);
			keywords_cells{kgc} = {DiStr,DStr,RStr};
			kgc = kgc+1;
		end
	end
end


%% ====================
% 2. Collect csv files and combine all data into a single structure 
% {'subnuclei','soma','D1','D2','D3','7','14','21'}
% keywords_region{1} = {'subnuclei'}; % The first element in the cell will be used to mark the combined data
% keywords_region{2} = {'soma'}; % The first element in the cell will be used to mark the combined data

% Group csv data according to keywords below

csv_folder = uigetdir(DefaultDir.csv, 'Select a folder to read csv files in it');
if csv_folder ~= 0
	DefaultDir.csv = csv_folder;
end

% combined_csv_data = cell(size(keywords_region));
% RegionType_num = numel(keywords_region);
for gn = 1:group_num
	keywords = keywords_cells{gn};

	[combined_csv_data{gn},csv_num,csv_list] = read_csv_files_in_folder(csv_folder,'keywords',keywords,...
		'gui_read',false);
	if csv_num ~= 0
		combined_csv_data{gn}.label = keywords_groups{gn};
		combined_csv_data{gn}.dilution = keywords{1};
		combined_csv_data{gn}.days = keywords{2};
		combined_csv_data{gn}.region = keywords{3};

		% keywords_WithSpace = cellfun(@(x) [x, ' '], keywords, 'UniformOutput',false); % add space to the end of each keyword
		% keywords_string = string(keywords_WithSpace); % convert keywords cell array to string array for display
		% combined_csv_data{gn}.RegionType = sprintf('%s',keywords_string); % add keywords as label
	end
end 
roi_val_data = [combined_csv_data{:}]; % Put csv_data containing different region-type data into a single struct  

%% ====================
% 3. Save the roi_data
default_mat_path = fullfile(DefaultDir.mat,'*.mat');
mat_folder = uigetdir(DefaultDir.mat,'Select a folder to save fluorescence intensity measurement')
% [mat_file,mat_folder] = uiputfile(default_mat_path,'Save fluorescence intensity measurement', 'FIM.mat');
if mat_folder ~= 0
	DefaultDir.mat = mat_folder;
	dt = datestr(now, 'yyyymmdd-HHMM');
	save(fullfile(mat_folder, [dt, '_FIM_data']),...
	    'roi_val_data');
end

%% ====================
% 4. Organize data and calculate the mean, std and ste for each group 
grouped_data = empty_content_struct({'region','data'},R_num); % group data using region type ('subnuclei' or 'soma')
for rn = 1:R_num
	region = RegionGroup{rn};
	grouped_data(rn).region = region;
	[~,idx] = filter_CharCells({roi_val_data.region},region);
	region_data = roi_val_data(idx); % data belongs to certain region type: subnuclei or soma

	grouped_data(rn).data = empty_content_struct({'label','dilutionID','days','FijiTbl','RawVal','MeanVal','StdVal','SteVal'},...
		group_num_MonoRegion);

	for gnmr = 1:group_num_MonoRegion
		grouped_data(rn).data(gnmr).label = keywords_groups_MonoRegion{gnmr};
		label = keywords_groups_MonoRegion{gnmr};
		[~,idx] = filter_CharCells({region_data.label},label);
		grouped_data(rn).data(gnmr).FijiTbl = vertcat(region_data(idx).combined_data);
		FijiTbl = vertcat(region_data(idx).combined_data);

		if ~isempty(grouped_data(rn).data(gnmr).FijiTbl)
			grouped_data(rn).data(gnmr).dilutionID = region_data(idx(1)).dilution;
			grouped_data(rn).data(gnmr).days = region_data(idx(1)).days;

			grouped_data(rn).data(gnmr).RawVal = [FijiTbl.Mean];
			grouped_data(rn).data(gnmr).MeanVal = mean([FijiTbl.Mean]);
			grouped_data(rn).data(gnmr).StdVal = std([FijiTbl.Mean]);
			grouped_data(rn).data(gnmr).SteVal = std([FijiTbl.Mean]) / sqrt(length([FijiTbl.Mean]));
		end
	end 
end

%% ====================
% 5. Plot data
close all
FontSize = 18;
FontWeight = 'bold';
plot_region = 'subnuclei'; % 'subnuclei' or 'soma'

% Get data measured using ROI type specified by 'plot_region'
[~,idx_region] = filter_CharCells({grouped_data.region},plot_region);
region_data = grouped_data(idx_region).data;

% Group data according to the dilution and make box plots showing the results of various waiting days
Line_y = cell(1,Di_num);
Line_y_error = cell(1,Di_num);
Line_xlabel = cell(1,Di_num);
[f_box] = fig_canvas(Di_num,'fig_name','BoxPlot SameDilution DifferentWaitingTime'); % create figure and customize the size
tlo = tiledlayout(f_box, ceil(Di_num/3), 3); % set tiles for plots
for din = 1:Di_num
	ax = nexttile(tlo);
	dilutionID = DilutionGroup{din};

	[~,idx_di] = filter_CharCells({region_data.dilutionID},dilutionID);
	di_data = region_data(idx_di);
	Line_y{din} = [di_data.MeanVal]; % collect data for line plot
	Line_y_error{din} = [di_data.SteVal]; % collect data for line plot
	% Line_xlabel{din} = {di_data.days}; % collect data for line plot

	RawData_di = {di_data.RawVal};
	groupNames = {di_data.label};
	[~, box_stat.(din)] = boxPlot_with_scatter(RawData_di, 'groupNames', groupNames,...
		'plotWhere', ax, 'stat', true, 'FontSize', FontSize,'FontWeight',FontWeight);
end

% Use line plot to show the fluorecense levels at different time. One line one dilution. All dilutions are plotted in one axis
[f_line] = fig_canvas(1,'fig_name','LinePlot');
ax = gca;
Line_x = [1:1:D_num];
Line_xlabel = DayGroup;
hold on
for din = 1:Di_num
	errorbar(Line_x,Line_y{din},Line_y_error(din));
end
xticks(Line_x);
xticklabels(Line_xlabel);