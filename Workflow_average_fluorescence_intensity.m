% This workflow can load csv files containing FIJI measurement in chosen folder, group the data 
% according to dilution, waiting-time, and ROI type, and finally plot the results  
%
% Pre-requirements (find details in "FIJI workflow for fluorescence intensity measurement" on Overleaf)
% 	- Correct background 
% 	- Correct shading (optional): no need to do this if the illumination is even 
% 	- Use FIJI/ImageJ to measure fluorescence intensity of ROIs in corrected image and save the measurements in csv files
% 
% 
% First of all: If you are using this code for the first time, modify some default settings below
% - Modify the FIM_folder (fluorescence intensity measurement).
%	FIM_folder will be the default location for loading csv files and save data. You will start to navigate from here in GUI
PC_name = getenv('COMPUTERNAME'); 
% Change the PC_name in the following lines, such as 'GD-AW-OFFICE' or 'LAPTOP-84IERS3H' to your own computer's name, and specify your FIM_folder
if strcmp(PC_name, 'GD-AW-OFFICE') % GD's office desktop
	FIM_folder = 'F:\Workspace\Confocal'; 
elseif strcmp(PC_name, 'LAPTOP-84IERS3H') % GD's laptop
	FIM_folder = 'C:\Users\guoda\Documents\Workspace\Analysis\Confocal'; 
end
%
%
% - Modify keywords: Keywords are used in csv file names. They represent dilution group, waiting time, and ROI types
% Example csv file name: ****_D1_07_soma-1_**** (dilution1, 7 waiting days, soma ROIs, first csv file in this group)
DilutionGroup = {'D1','D2','D3'}; % Dilution1, Dilution2, etc. 
DayGroup = [7 14 21]; % Waiting time in days
RegionGroup = {'subnuclei','soma'}; % ROI type


%% ====================
% 1. Prepare to analysis

% Assign default folders
DefaultDir.csv = FIM_folder;
DefaultDir.mat = FIM_folder;
DefaultDir.fig = FIM_folder;

% Get the number of keywords in each category
Di_num = numel(DilutionGroup);
D_num = numel(DayGroup);
R_num = numel(RegionGroup);

% Group keywords 
group_num = Di_num*D_num*R_num;
keywords_groups = cell(1,group_num); % each cell contains a single string, such as 'D1-07-soma'
% keywords_cells = cell(1,group_num); % each cell contains a cell of keywords, such as {'D1','14','soma'}
combined_csv_data = cell(1,group_num);
kgc = 1; % current keyword group counting
group_num_MonoRegion = Di_num*D_num;
keywords_groups_MonoRegion = cell(1,group_num_MonoRegion);
kgc_mr = 1;
for din = 1:Di_num
	DiStr = DilutionGroup{din};
	for dn = 1:D_num
		Dnum = DayGroup(dn);
		DStr = num2str(Dnum,'%02.f');
		keywords_groups_MonoRegion{kgc_mr} = sprintf('%s-%s', DiStr,DStr); % example: 'D1-07'
		for rn = 1:R_num
			RStr = RegionGroup{rn}; % example: 'soma' or 'nuclei'
			keywords_groups{kgc} = sprintf('%s-%s-%s',DiStr,DStr,RStr);
			keywords_cells_1{kgc} = {DiStr,DStr,RStr};
			keywords_cells_2{kgc} = {keywords_groups_MonoRegion{kgc_mr},RStr};
			kgc = kgc+1;
		end
		kgc_mr = kgc_mr+1;
	end
end


%% ====================
% 2. Collect csv files and combine all data into a single structure 

% Select the folder to load csv files with GUI
csv_folder = uigetdir(DefaultDir.csv, 'Select a folder to read csv files in it');
if csv_folder ~= 0
	DefaultDir.csv = csv_folder;
end

% Load csv files and save the data in a struct var
for gn = 1:group_num
    keywords_1 = keywords_cells_1{gn};
    keywords_2 = keywords_cells_2{gn};

	[combined_csv_data{gn},csv_num,csv_list] = read_csv_files_in_folder(csv_folder,'keywords',keywords_2,...
		'del_col',{'StdDev'},'gui_read',false,'debug_mode',true);
	if csv_num ~= 0
		combined_csv_data{gn}.label = keywords_groups{gn};
		combined_csv_data{gn}.dilution = keywords_1{1};
		combined_csv_data{gn}.days = str2double(keywords_1{2});
		combined_csv_data{gn}.region = keywords_1{3};
	end
end 
roi_val_data = [combined_csv_data{:}]; % Put csv_data containing different region-type data into a single struct  

%% ====================
% 3. Save the roi_val_data with GUI
default_mat_path = fullfile(DefaultDir.mat,'*.mat');
mat_folder = uigetdir(DefaultDir.mat,'Select a folder to save fluorescence intensity measurement');
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
	[~,idx_region] = filter_CharCells({roi_val_data.region},region);
	region_data = roi_val_data(idx_region); % data belongs to certain region type: subnuclei or soma

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
		else
			grouped_data(rn).data(gnmr).dilutionID = '';
			grouped_data(rn).data(gnmr).days = [];
		end
	end 
end

%% ====================
% 5. Plot data
close all
plot_region = 'soma'; % 'subnuclei' or 'soma'. Specify the ROI type

% Settings for figures
SaveFig = false; % true/false
FontSize = 18;
FontWeight = 'bold';
LineColors = {'#3FF5E6', '#F55E58', '#F5A427', '#4CA9F5', '#33F577',...
		'#408F87', '#8F4F7A', '#798F7D', '#8F7832', '#28398F', '#000000'};
LineWidth = 1;

% Get data measured using ROI type specified by 'plot_region'
[~,idx_region] = filter_CharCells({grouped_data.region},plot_region);
region_data = grouped_data(idx_region).data;

% Group data according to the dilution and make box plots showing the results of various waiting days
Line_x = cell(1,Di_num);
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
    Line_x{din} = [di_data.days];
	Line_y{din} = [di_data.MeanVal]; % collect data for line plot
	Line_y_error{din} = [di_data.SteVal]; % collect data for line plot
	% Line_xlabel{din} = {di_data.days}; % collect data for line plot

	RawData_di = {di_data.RawVal};
	groupNames = {di_data.label};
	[~, box_stat.(DilutionGroup{din})] = boxPlot_with_scatter(RawData_di, 'groupNames', groupNames,...
		'plotWhere', ax, 'stat', true, 'FontSize', FontSize,'FontWeight',FontWeight);
end
LegendStr = DilutionGroup; % add dilution name as legend for line plots

% Use line plot to show the fluorecense levels at different time. One line one dilution. All dilutions are plotted in one axis
[f_line] = fig_canvas(1,'fig_name','LinePlot of Multiple Dilutions');
ax = gca;
% Line_x = [1:1:D_num];
Line_xlabel = DayGroup;
hold on
for din = 1:Di_num
	errorbar(Line_x{din},Line_y{din},Line_y_error{din},...
		'Color',LineColors{din},'LineWidth',LineWidth);
end
xlim([0 28])
xticks(DayGroup);
xticklabels(arrayfun(@num2str, DayGroup, 'UniformOutput', 0));
legend(LegendStr,'location','northeast')
legend('boxoff')
set(gca,'box','off')
set(gca,'FontSize',FontSize)
set(gca,'FontWeight',FontWeight)


if SaveFig
	% default_fig_path = fullfile(DefaultDir.fig,'*.mat');
	FigFolder = uigetdir(DefaultDir.fig,'Select a folder to save figures');
	% [mat_file,mat_folder] = uiputfile(default_mat_path,'Save fluorescence intensity measurement', 'FIM.mat');
	if mat_folder ~= 0
		DefaultDir.fig = FigFolder;
		dt = datestr(now, 'yyyymmdd-HHMM');
		fbox_name = sprintf('%s_%s_boxplot',dt,plot_region);
		savePlot(f_box,...
			'guiSave', 'off', 'save_dir', FigFolder, 'fname', fbox_name);

		fline_name = sprintf('%s_%s_lineplot',dt,plot_region);
		savePlot(f_line,...
			'guiSave', 'off', 'save_dir', FigFolder, 'fname', fline_name);

		save(fullfile(FigFolder, [dt, '_boxplot_stat']),...
		    'box_stat');
	end
end