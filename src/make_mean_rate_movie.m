% Script to generate viterbi movie illustrating mean mRNA dynamics
% should be run on the server
%%% Load Data & Create Write Paths
clear 
close all
%Embryo ID and Folder location
Prefix = '2018-07-19-Eve2Ms2KruppelLlamaA_01';
PrefixDropboxFolder = '../dat/Kruppel_eve2_pass1';
FISHPath = '../../figure_data/2017-06-15-eve2_20sec_10_raw/';
% WritePath = '../../figures/Fig4/viterbi_movie_frames/';
% mkdir(WritePath);
%Load the data
load([PrefixDropboxFolder,'\CompiledParticles.mat'])
load([PrefixDropboxFolder,'\' Prefix,'_lin.mat'])
load([PrefixDropboxFolder,'\Ellipses.mat'])

%%% ------ movie parameters
MaxRadius = 22;
[px, py] = meshgrid(1:512,1:256);
MaxmRNAFluo = 1000;
CurrentNC = 14;
p_fluo = 12.74; % Specify AU - PolII Calibration
fill_color = [85 169 116]/256;%[122 169 116]/256; % RGB for patch fills
elongation_time = (140*(1-1302/6544) + .5*140*1302/6544)/60; % minutes 
%%% ------- ID Variables
project = 'mHMMeve2_weka_inf_2018_05_07'; %project identifier
%%% write path
WritePath = '..\..\figures\Fig3\';
mkdir([WritePath '/MeanRateFrames/'])
% load trace data
load(['../../figure_data/inference_traces_' project '_dT20.mat'])
%%% ------------------- Generate Movies ------------------------------- %%%
setID = 10;
last_frame = max([trace_struct_final([trace_struct_final.setID]==setID).all_frames]);
set_struct = trace_struct_final([trace_struct_final.setID]==setID);
% first_frame = min([trace_struct_filtered([trace_struct_filtered.setID]==setID).all_frames]);
FrameRange=nc14:last_frame;

n_ticks = 5; % set tick increment
%Iterate Through Frames
for Dummy = FrameRange    
    CurrentFrame = max(FrameRange) - Dummy + nc14;
    %Highlight Active Regions with Viterbi State    
    %Track pixel assignments
    NucleusStateMat = zeros(256,512,3);
    NucleusIDMat = zeros(256,512);
    NucleusDistMat = ones(256,512)*MaxRadius;
    %Loop through ALL nuclei (inactive and active) and assign patches
    for s = 1:length(schnitzcells)
        MaxFrame = max(schnitzcells(s).frames);
        MinFrame = min(schnitzcells(s).frames);        
        if sum(ismember(MaxFrame,120:175)) == 1 && CurrentFrame > MaxFrame
            cf = MaxFrame;
        else
            cf = CurrentFrame;
        end
        CurrentEllipse= schnitzcells(s).cellno(...
                        schnitzcells(s).frames==...
                        cf);  
        x =  Ellipses{cf}(CurrentEllipse,1)+1;
        y =  Ellipses{cf}(CurrentEllipse,2)+1;
        % Exclude Nuclei that start midway through
        if sum(ismember(MinFrame,120:175)) 
            continue
        end
        if isempty(x)
            continue
        end 
        distances = ((px-x).^2 + (py-y).^2).^.5;
        candidate_indices = NucleusDistMat > distances; 
        %Record Fluorescence        
        NucleusIDMat(candidate_indices) = s;
        NucleusDistMat(candidate_indices) = distances(candidate_indices);
    end
    %Loop Through Particles to see which remain on and how long these have
    %been active
    ParticlesToShow = [];
    for i = 1:length(set_struct)        
        all_frames = set_struct(i).all_frames;                     
        cp_frames = set_struct(i).cp_frames;                     
        extant_frames = all_frames((all_frames <= CurrentFrame)&...
            (all_frames>nc14));        
        if  CurrentFrame < nc14
            continue        
        else
            if ~isempty(extant_frames)               
%                 try     
                    %Record Fluorescence                    
                    filter = NucleusIDMat==set_struct(i).Nucleus;
                    frame_filter = ismember(cp_frames,extant_frames);
                    mf = nanmean(set_struct(i).fluo(frame_filter));
                    mf = min(MaxmRNAFluo,mf)/MaxmRNAFluo;
                    for k = 1:3
                        slice = NucleusStateMat(:,:,k);                        
                        slice(filter) = fill_color(k)*mf;
                        NucleusStateMat(:,:,k) = slice;
                    end                    
                    ParticlesToShow = [ParticlesToShow set_struct(i).Nucleus];
%                 catch
%                     display(['Error in frame ',num2str(CurrentFrame)])                
%                 end
            end
        end
    end
    %Now Draw Nucleus Borders for active nuclei
    NucleusBorderMat = zeros(size(NucleusIDMat));
    window = 1;
%     c_filter = ones(5,5);
    for i = ParticlesToShow
        %get coordinates of nucleus patch
        if sum(sum(NucleusIDMat==i)) > 0
            x_vec = reshape(px(NucleusIDMat==i),[],1);
            y_vec = reshape(py(NucleusIDMat==i),[],1);
            for j = 1:length(x_vec)
                metric = sum(sum(NucleusIDMat(max(1,y_vec(j)-window):min(y_vec(j)+window,256),...
                         max(1,x_vec(j)-window):min(x_vec(j) + window,512))));
                if metric~= i*(2*window+1)^2
                    NucleusBorderMat(y_vec(j),x_vec(j)) = 1;
                end
            end
        end        
    end
    %Prevent overlap between fluroescence mask and borders
    for k = 1:3
        slice = NucleusStateMat(:,:,k);
        slice(NucleusBorderMat>0) = 0;
        NucleusStateMat(:,:,k) = slice;
    end

    %Make a maximum projection of the mRNA channel
    D=dir([FISHPath,filesep,Prefix,'_',num2str(CurrentFrame,'%03d'),'_z*.tif']);
    %Do not load the first and last frame as they are black
    ImageTemp=[];
    for i=2:(length(D)-1)
        ImageTemp(:,:,i-1)=imread([FISHPath,filesep,D(i).name]);
    end
    mRNAImageRaw=max(ImageTemp,[],3);    
%         mRNAImage = mRNAImage/max(mRNAImage(:));
    %Load the corresponding histone image
    HistoneImage=imread([FISHPath,filesep,...
        Prefix,'-His_',num2str(CurrentFrame,'%03d'),'.tif']);        

    %Overlay all channels
%     NucleusStateMat(NucleusStateMat>MaxmRNAFluo) = MaxmRNAFluo; % enforce 1000AU maximum
%     MCPshading = NucleusStateMat/MaxmRNAFluo;
%     MCPshading(MCPshading>1) = 1;
    MCPshading = NucleusStateMat(:,:,2);
    MCPChannel =  mRNAImageRaw/150 + MCPshading; %
    MCPChannel(MCPChannel>1) = 1;
    
    if CurrentFrame == max(FrameRange)
        MCPColorbar = MCPshading;
    end        
    HistoneChannel=  mat2gray(HistoneImage) + NucleusStateMat(:,:,1);%
    HistoneChannel(HistoneChannel>1) = 1;
    StateChannel = NucleusStateMat(:,:,3);
    if CurrentFrame == max(FrameRange)
        ColorbarOverlay = NucleusStateMat;
        max_fluo = MaxmRNAFluo/p_fluo/elongation_time;
        tick_increment = round(max_fluo/n_ticks/3,1);
        tick_string = {'0'};
        for i = 1:n_ticks
            tick_string = [tick_string{:} {num2str(round(i*tick_increment))}];
        end
        cb_struct_mr.ColorbarOverlay = ColorbarOverlay;
        cb_struct_mr.n_ticks = n_ticks;
        cb_struct_mr.tick_string = tick_string;
        cb_struct_mr.tick_increment = tick_increment;
        cb_struct_mr.max_fluo = max_fluo;
        save([WritePath 'MeanRateColorBarInfo.mat'],'cb_struct_mr')                                
    end
    ImOverlay=cat(3,HistoneChannel,MCPChannel,StateChannel);
    % Make Vectors for Colormap    
    OverlayFig = figure(1);%figure('Visible','off');
    clf   
 
    imshow(fliplr(ImOverlay),'DisplayRange',[])
    
    ylim([16,240])
    if ceil(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))) == 40        
        saveas(gcf,[WritePath '\mean_rate_patch_overlay_40min.tif']); 
    end
%     h = colorbar('XTick', (0:(100/n_ticks):100)/100, 'XTickLabel', tick_string);
%     h.Label.String = 'activity duration (minutes)';
    text(446,225,[num2str(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),'%02d'),...
        ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])
%     xlim([0,512])    
    
    drawnow
    saveas(gcf,[WritePath '\MeanRateFrames\nc',num2str(CurrentNC),...
        '-',num2str(CurrentFrame,'%02d'),'.tif']);       
end
