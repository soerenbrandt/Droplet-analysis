classdef videoMaker
    % videoMaker V1.52
    % added adjustments for number average and volume average to
    % plotDistribution (not yet fully implemented)
    % simplified volume calculation
    % fixed bug that occured when obj.setupRegion was called with an empty
    % obj.rect
    
    % Analyses objects in video and returns data or video with
    % respective features and analysis.
    %
    % Object detection (based on videoFrame class) includes
    %   'circle'                     - detects circular objects in video frame
    %   'bright' (default) or 'dark'  - object polarity compared to background
    %   'fluorescence' (default) or 'brightfield' - clarifies image mode
    %   'rect'                        - specify image region for detection
    %
    % Analysis options include
    %   'tracking'      - returns particle tracks from image
    %   'streamlines'   - calculates local flow velocities (includes tracking) 
    %   'distribution'  - returns particle distribution, include regions to
    %                     specify regions of interest
    %   'regions'       - specify regions of interest
    %   Analysis options can be specified during construction of videoMaker
    %   class or during makeVideo to include data in stream
    
    properties
        frames % array with analysis of each frame
        videoLink % file destination to original video
        analysis % = struct('tracking',[],'distribution',[],'streamlines',[]); % analysis from frame-by-frame data
        imageType % image appearance ('fluorescence', 'brightfield', or 'darkfield')
    end
    properties (Dependent)
        Objects % objects found in image analysis
        Distribution % distribution of particles in video
    end
    properties (SetAccess = protected)
        forType % analysis performed for object types
        forAnalysis % analysis methods performed
        rect % rectangle to crop the video
        regions % relevant regions in image
        options % additional parameters set for analysis and image processing
        scale % for converting pixel units to physical units (factor in um/pixel, depth in um)
    end
    properties (SetAccess = immutable)
        date % creation date of the video file
    end
    
    %% Methods
    methods
        % Video analysis functions
        function obj = videoMaker(videoLink,varargin) %WRITE
            % obj = videoMaker(videoLink,varargin) creates a videoMaker
            % object from tracking objects in the videoFile specified by 
            % videoLink. varargin can be used to define a custom sequence of
            % analysis tools, including 'circle', 'tracking', and 'Show:
            % <per number of frames>'. Default is varargin =
            % {'rescale','circle','tracking','show: 500'}.
            %
            %   obj = videoMaker() calls the videoMaker class with the
            %   standard sequence {'rescale','circle','tracking','show: 500
            %   '}.
            %
            %   See also videoMaker.runAnalysis, videoMaker.trackParticles,
            %   and videoFrame.
            
            %---------------------------------------------------------------
            %         Step 1 Set default variables
            %---------------------------------------------------------------
            obj.date = date;
            if nargin < 1 || isempty(videoLink)
                fileTypes = {'*.m4v','*.mov','*.avi','*.mp4'};
                [name, path] = uigetfile({strjoin(fileTypes,';'),'All Video Files'},'','~/Desktop/');
                if ~ischar(name)
                    return
                end
                obj.videoLink = [path, name];
            else
                obj.videoLink = videoLink;
            end
            
            if any(cellfun(@(str)strcmp(str,'odyssee'),lower(varargin)))
                obj = obj.setupVideoAnalysisFor(varargin{:}, 'full');
                obj.prepareForOdyssee;
                return
            elseif isempty(varargin)
                varargin = {'rescale','circle','tracking','show: 500'}; % standard operations
            elseif ~any(ismember({'circle','full'},lower(varargin)))
                varargin{end+1} = 'circle';
            end
            
            %---------------------------------------------------------------
            %         Step 2 Set analysis variables
            %---------------------------------------------------------------
            obj = obj.setupVideoAnalysisFor(varargin{:});
            
            %---------------------------------------------------------------
            %         Step 3 Analyse frames
            %---------------------------------------------------------------
            obj.frames = obj.analyseFrames;
            
            %---------------------------------------------------------------
            %         Step 4 Analyse video
            %---------------------------------------------------------------
            obj.analysis = obj.analyseVideo;
        end
        function obj = redoAnalysis(obj,varargin)
            % obj = redoAnalysis(obj,varargin) repeats the analysis for obj
            % or performs the analysis specified in varargin.
            %
            %   obj = redoAnalysis(obj,'tracking') performs particle
            %   tracking analaysis for obj.
            %
            %   See also videoMaker.runAnalysis, videoMaker.trackParticles.
            
            if isempty(varargin)
                varargin = [obj.forType, obj.forAnalysis, obj.options, {obj.imageType}];
            end
            
            %---------------------------------------------------------------
            %         Step 1 Set analysis variables
            %---------------------------------------------------------------
            obj = setupVideoAnalysisFor(obj,varargin{:});
            
            %---------------------------------------------------------------
            %         Step 2 Analyse frames
            %---------------------------------------------------------------
            if ~isempty(obj.forType)
                obj.frames = analyseFrames(obj);
            end
            
            %---------------------------------------------------------------
            %         Step 3 Analyse video
            %---------------------------------------------------------------
            obj.analysis = analyseVideo(obj);
        end % check which analysis has been done and redo (make different from runAnalysis
        function analysis = runAnalysis(obj,of)
            % analysis = runAnalysis(obj,of) returns the analysis specified
            % by of for the videoMaker object obj.
            %
            %   analysis = runAnalysis(obj,'tracking') performs particle
            %   tracking analaysis for obj.
            %
            %   See also videoMaker.runAnalysis, videoMaker.trackParticles.
            
            switch of
                case 'tracking'
                    analysis = trackParticles(obj);
                otherwise
                    analysis = [];
                    return
            end
        end
        
        % Storage and output functions
        function makeVideo(obj,varargin) % UPDATE, quick makeVideo, and pretty makeVideo + 
            % HELP for makeVideo currently unavailable.
            
            progress = waitbar(0,'checking analysis');
            %---------------------------------------------------------------
            %         Step 1 Check for requested elements
            %---------------------------------------------------------------
%             plots = {};
%             if ismember(varargin,'tracking') && isfield(obj.analysis,'tracking')
%                 plots = {@(obj,ax,Nr)videoMaker.plotTracks(obj,ax,Nr)};
%             end
            
            %---------------------------------------------------------------
            %         Step 2 Make video
            %---------------------------------------------------------------
            waitbar(0,progress,'locating video');
            video = obj.openVideo;
            
            % make video frame
            fig = figure; title('video frame');
            ax = axes(fig);
            firstFrame = imcrop(readFrame(video),obj.rect); imshow(firstFrame);
            ax.XLimMode = 'manual'; ax.XTickMode = 'manual'; ax.YLimMode = 'manual'; ax.YTickMode = 'manual';
%             function [ax,fig] = makeFrame(video,obj,frameInd)
%                 time = video.CurrentTime;
%                 fig = figure;
%                 image(imcrop(readFrame(video),obj.rect)); ax = fig.Children;
%                 set(fig, 'Position', [0 0 obj.rect(3:4)*1.0526]);
%                 set(fig.Children, 'Position', [0.025 0.025 0.95 0.95],'Units','normalized');
%                 ax.XLimMode = 'manual'; ax.XTickMode = 'manual'; ax.YLimMode = 'manual'; ax.YTickMode = 'manual';
%                 imshow(readFrame(video));
%                 ax = fig.Children;
%                 if time ~= obj.frames(frameInd).time
%                     box(ax,'on'); ax.XTick = []; ax.YTick = []; 
%                     ax.XColor = 'r'; ax.YColor = 'r'; ax.LineWidth = 2; 
%                 end
%             end
            
            framesTotal = length(obj.frames); timeRem = nan(1,framesTotal);
            frameInd = 2; % editedFrames(numel(obj.frames)) = getframe(makeFrame(video,obj,1)); editedFrames = fliplr(editedFrames); video.CurrentTime = 0;
            editedFrames(numel(obj.frames)) = im2frame(firstFrame); editedFrames = fliplr(editedFrames); % video.CurrentTime = 0;
            while hasFrame(video) && video.CurrentTime < video.Duration*1.05
                tic % start clock
                time = video.CurrentTime;
                waitbar(time/video.Duration,progress,...
                    ['writing frame ',num2str(frameInd),' of ',num2str(framesTotal),...
                    ' (time rem: ',num2str(round((framesTotal-frameInd)*nanmean(timeRem)/60)),'min)']);
                
                % [ax,fig] = makeFrame(video,obj,frameInd);
                imshow(imcrop(readFrame(video),obj.rect)); % was image before, faster?
                ax = fig.Children;
                
                % PLOT DATA FROM CELL ARRAY OF FUNCTION HANDLES
%                 for add = plots
%                     add{1}(obj,ax,frameInd); % TAKES TOO LONG
%                 end
                
%                 if time ~= obj.frames(frameInd).time
%                     box(ax,'on'); ax.XTick = []; ax.YTick = []; 
%                     ax.XColor = 'r'; ax.YColor = 'r'; ax.LineWidth = 2; 
%                 end
                videoMaker.plotTracks(obj,ax,frameInd);
%                 currentFrame = videoMaker.addTracksToFrame(obj,frameInd,imcrop(readFrame(video),obj.rect),[1 1 1]);

                editedFrames(frameInd) = getframe(ax); %im2frame(currentFrame); %getframe(ax);
                % cla(ax); close(fig);
                frameInd = frameInd + 1;
                timeRem(frameInd) = toc;
            end
%             close(fig)

            waitbar(1,progress,'writing video file');
            
            [path,name] = fileparts(obj.videoLink);
            
            newVideo = VideoWriter([path,'/',name,'-analyzed'],'MPEG-4'); % CHANGE in varargin?
            newVideo.FrameRate = video.FrameRate; % CHANGE in varargin?

            % Clean out empty editedFrames
            editedFrames(arrayfun(@(x)isempty(x.cdata),editedFrames)) = [];

            % Write video file
            open(newVideo)
            writeVideo(newVideo,editedFrames);
            close(newVideo)
            
            delete(progress)
        end % FINISH
        function makeCrudeVideo(obj,varargin) % UPDATE, quick makeVideo, and pretty makeVideo + 
            % HELP for makeCrudeVideo currently unavailable.
            
            progress = waitbar(0,'checking analysis');
            %---------------------------------------------------------------
            %         Step 1 Check for requested elements
            %---------------------------------------------------------------
            plots = {};
            if ismember(varargin,'tracking') && isfield(obj.analysis,'tracking')
                plots = {@(obj,ax,Nr)videoMaker.plotTracks(obj,ax,Nr)};
            end
            
            %---------------------------------------------------------------
            %         Step 2 Make video
            %---------------------------------------------------------------
            waitbar(0,progress,'locating video');
            video = obj.openVideo;
            
            % function [ax,fig] = makeFrame(video,obj,frameInd)
                %time = video.CurrentTime;
%                 fig = figure; set(fig, 'Visible', 'off');
%                 image(imcrop(readFrame(video),obj.rect)); ax = fig.Children;
%                 set(fig, 'Position', [0 0 obj.rect(3:4)*1.0526]);
%                 set(fig.Children, 'Position', [0.025 0.025 0.95 0.95],'Units','normalized');
%                 ax.XLimMode = 'manual'; ax.XTickMode = 'manual'; ax.YLimMode = 'manual'; ax.YTickMode = 'manual';
%                 imshow(readFrame(video));
%                 ax = fig.Children;
%                 if time ~= obj.frames(frameInd).time
%                     box(ax,'on'); ax.XTick = []; ax.YTick = []; 
%                     ax.XColor = 'r'; ax.YColor = 'r'; ax.LineWidth = 2; 
%                 end
            % end
            
            framesTotal = length(obj.frames); timeRem = nan(1,framesTotal);
            frameInd = 1; % editedFrames(numel(obj.frames)) = getframe(makeFrame(video,obj,1)); editedFrames = fliplr(editedFrames); video.CurrentTime = 0;
            editedFrames(numel(obj.frames)) = im2frame(readFrame(video)); editedFrames = fliplr(editedFrames); video.CurrentTime = 0;
            while hasFrame(video) && video.CurrentTime < video.Duration*1.05
                tic % start clock
                time = video.CurrentTime;
                waitbar(time/video.Duration,progress,...
                    ['writing frame ',num2str(frameInd),' of ',num2str(framesTotal),...
                    ' (time rem: ',num2str(round((framesTotal-frameInd)*nanmean(timeRem)/60)),'min)']);
                
                % [ax,fig] = makeFrame(video,obj,frameInd);
                %image(imcrop(readFrame(video),obj.rect));
                
                % PLOT DATA FROM CELL ARRAY OF FUNCTION HANDLES
%                 for add = plots
%                     add{1}(obj,ax,frameInd); % TAKES TOO LONG
%                 end
                
%                 if time ~= obj.frames(frameInd).time
%                     box(ax,'on'); ax.XTick = []; ax.YTick = []; 
%                     ax.XColor = 'r'; ax.YColor = 'r'; ax.LineWidth = 2; 
%                 end
                currentFrame = videoMaker.addTracksToFrame(obj,frameInd,imcrop(readFrame(video),obj.rect),[1 1 1]);

                editedFrames(frameInd) = im2frame(currentFrame); %getframe(ax);
%                 cla(ax); %close(fig);
                frameInd = frameInd + 1;
                timeRem(frameInd) = toc;
            end
%             close(fig)

            waitbar(1,progress,'writing video file');
            
            [path,name] = fileparts(obj.videoLink);
            
            newVideo = VideoWriter([path,'/',name,'-analyzed'],'MPEG-4'); % CHANGE in varargin?
            newVideo.FrameRate = video.FrameRate; % CHANGE in varargin?

            % Clean out empty editedFrames
            editedFrames(arrayfun(@(x)isempty(x.cdata),editedFrames)) = [];

            % Write video file
            open(newVideo)
            writeVideo(newVideo,editedFrames);
            close(newVideo)
            
            delete(progress)
        end % FINISH
        function trackObject(obj,n)
            % trackObject(obj,n) returns an image sequence of the video
            % analyzed by obj. The image sequence tracks one particle if n
            % is an integer, or all particles in frames specified by n if n
            % is an array of integers.
            %
            %   See also videoMaker.trackParticles.
            
            if nargin < 2
                return
            end
            
            % open video
            video = obj.openVideo;
            
            switch length(n)
                case 1
                    try object = obj.Objects(n);
                    catch me
                        error(me)
                    end
                    
                    startTime = max([object.atTime(1)-0.4,0]); endTime = object.atTime(end)+0.4;
                    
                    % select relevant frames and highlight object
                    video.CurrentTime = startTime;
                    
                    numberOfFrames = round((endTime-startTime)*video.FrameRate);
                    sequence = uint8(zeros(video.Height,video.Width,3,numberOfFrames));
                    
                    frame = 1;
                    while hasFrame(video) && video.CurrentTime <= endTime
                        location = find(video.CurrentTime == object.atTime);
                        Im = readFrame(video);
                        sequence(:,:,:,frame) = insertShape(Im,'circle', ...
                            [object.position(location,:)+obj.rect(1:2) object.size(location)/obj.scale.factor], ...
                            'LineWidth',2,'Color','red');
                        frame = frame + 1;
                    end
                otherwise
                    % select relevant frames and highlight objects
                    objects = obj.Objects(arrayfun(@(object)any(ismember(n,object.inFrame)),obj.Objects));
                    colors = uint8(255*hsv(length(objects)));
                    
                    % remove n larger than number of frames
                    if any(n > numel(obj.frames))
                        warning('Some indices excede number of frames');
                        n(n > numel(obj.frames)) = [];
                    end
                    
                    numberOfFrames = numel(n);
                    sequence = uint8(zeros(video.Height,video.Width,3,numberOfFrames));
                    
                    progress = waitbar(0,'grabbing frames');
                    for Nr = 1:numberOfFrames
                         waitbar(Nr/numel(n),progress,...
                                ['analyzing frame ',num2str(Nr),' of ',num2str(numberOfFrames)]);
                            
                        video.CurrentTime = obj.frames(n(Nr)).time;
                        locations = arrayfun(@(object)[object.position(video.CurrentTime == object.atTime,:) + obj.rect(1:2) ...
                                                       median(object.size)/obj.scale.factor],objects,'uni',0);
                        clean = cellfun(@(x)size(x,2)==3,locations);
                        Im = readFrame(video);
                        sequence(:,:,:,Nr) = insertShape(Im,'circle',vertcat(locations{clean}),'LineWidth',2,'Color',colors(clean,:));
                    end
                    delete(progress)
                    
%                     startTime = obj.frames(min(n)).time; endTime = obj.frames(max(n)).time;
%                     
%                     % select relevant frames and highlight object
%                     video.CurrentTime = startTime;
%                     
%                     while hasFrame(video) && video.CurrentTime <= endTime && frame <= max(n)
%                         locations = arrayfun(@(object)[object.position(video.CurrentTime == object.atTime,:) + obj.rect(1:2) ...
%                                                        median(object.size)/obj.scale.factor],objects,'uni',0);
%                         clean = cellfun(@(x)size(x,2)==3,locations);
%                         Im = readFrame(video);
%                         sequence(:,:,:,frame) = insertShape(Im,'circle',vertcat(locations{clean}),'LineWidth',2,'Color',colors(clean,:));
%                         frame = frame + 1;
%                     end
            end
            implay(sequence,video.FrameRate)
        end
        function makeCVS(obj,varargin) %WRITE
            %HELP for makeCVS currently unavailable
        end
        function save(obj,path)
            % save(obj,path) saves the videoMaker object obj as a matfile
            % to the path specified in path. The name of the videoMaker
            % object is the name of the video specified in obj.videoLink.
            %
            %   save(obj) saves the videoMaker object under the path
            %   specified by videoLink.
            %
            %   See also videoMaker.
            
            if nargin < 2 || isempty(path)
                [path,name] = fileparts(obj.videoLink);
            else
                [~,name] = fileparts(obj.videoLink);
            end
            
            save([path,'/',name,'.mat'],'obj')
        end
        
        % Data access functions
        function inRegion = objectsInRegion(obj,varargin)
            % inRegion = objectsInRegion(obj,varargin) returns all objects
            % specified as a list of integers in varargin. If no region 
            % from obj.regions is selected, it returns objects for all 
            % regions.
            %
            %   See also videoMaker.setupRegions, videoMaker.trackParticles.
            
            if isempty(obj.regions)
                obj = obj.setupRegions;
            end
            
            if isempty(varargin)
                regionIDs = 1:numel(obj.regions);
            else
                regionIDs = varargin{1};
            end
            
            objects = obj.Objects;
            
            % filter objects based on location of droplets
            inRegion = cell(size(regionIDs));
            for n = regionIDs
                if isempty(obj.rect)
                    regionXs = obj.regions{n}(1:end/2); regionYs = obj.regions{n}(1+end/2:end);
                else
                    regionXs = obj.regions{n}(1:end/2) - obj.rect(1); regionYs = obj.regions{n}(1+end/2:end) - obj.rect(2);
                end
                inRegion{n} = objects(cellfun(@(pos)any(inpolygon(pos(:,1),pos(:,2),regionXs,regionYs)),{objects.position}));
            end
        end
        function s = scatterDistribution(obj,ax,at,varargin)
            if isempty(obj.regions)
                obj = obj.setupRegions;
            end
            if nargin < 3 || isempty(at)
                at = xlim(ax); at = at(2);
            end
            if ismember('volume',varargin(cellfun(@(c)isa(c,'char'),varargin)))
                param = 'volume';
            else
                param = 'size';
            end
            
            if isempty(varargin(cellfun(@(c)isa(c,'numeric'),varargin)))
                regionIDs = 1:numel(obj.regions);
            else
                regionIDs = varargin{1};
            end
            
            if isempty(regionIDs)
                return
            else
                objects = obj.Objects;
            end
            
            if nargin < 2 || isempty(ax)
                ax = axes(figure);
                colormap(ax,flipud(bone))
            end
            
            % filter objects based on location of droplets
            hold(ax,'on')
            for n = regionIDs
                if isempty(obj.rect)
                    regionXs = obj.regions{n}(1:2:end-1); regionYs = obj.regions{n}(2:2:end);
                else
                    regionXs = obj.regions{n}(1:2:end-1) - obj.rect(1); regionYs = obj.regions{n}(2:2:end) - obj.rect(2);
                end
                inRegion = objects(cellfun(@(pos)any(inpolygon(pos(:,1),pos(:,2),regionXs,regionYs)),{objects.position}));
                
                meanParam = arrayfun(@(obj)mean(obj.(param)),inRegion(arrayfun(@(obj)length(obj.(param)),inRegion)>2));
                [values, edges] = histcounts(meanParam,100,'Normalization','Probability');
                s = scatter(ax,at*ones(size(values)),edges(2:end),[],values,'filled');
            end
            hold(ax,'off')
        end
        function h = plotHistogram(obj,ax,varargin)
            if isempty(obj.regions)
                obj = obj.setupRegions;
            end
            if ismember('volume',varargin(cellfun(@(c)isa(c,'char'),varargin)))
                param = 'volume';
            else
                param = 'size';
            end
            
            if isempty(varargin(cellfun(@(c)isa(c,'numeric'),varargin)))
                regionIDs = 1:numel(obj.regions);
            else
                regionIDs = varargin{1};
            end
            
            if isempty(regionIDs)
                return
            else
                objects = obj.Objects;
            end
            
            if nargin < 2 || isempty(ax)
                ax = axes(figure);
            end
            
            % filter objects based on location of droplets
            hold(ax,'on')
            colorOrder = get(ax, 'ColorOrder');
            for n = regionIDs
                if isempty(obj.rect)
                    regionXs = obj.regions{n}(1:2:end-1); regionYs = obj.regions{n}(2:2:end);
                else
                    regionXs = obj.regions{n}(1:2:end-1) - obj.rect(1); regionYs = obj.regions{n}(2:2:end) - obj.rect(2);
                end
                % inRegion = objects(cellfun(@(pos)any(inpolygon(pos(:,1),pos(:,2),regionXs,regionYs)),{objects.position}));
                inRegion = @(object)inpolygon(object.position(:,1),object.position(:,2),regionXs,regionYs);
                
%                 if n == 1 % added for AG1B-5
%                     inRegion(arrayfun(@(obj)mean(obj.size),inRegion)<50) = [];
%                 end
                meanParam = arrayfun(@(object)mean(object.(param)(inRegion(object))),objects(arrayfun(@(object)sum(inRegion(object))>0,objects)));
                %meanParam = cell2mat(arrayfun(@(object)object.(param)(inRegion(object)),objects(arrayfun(@(object)sum(inRegion(object))>2,objects)),'uni',0));
                %2 meanParam = cell2mat(arrayfun(@(obj)(obj.(param)),inRegion(arrayfun(@(obj)length(obj.(param)),inRegion)>=1),'uni',0));
                % meanParam = arrayfun(@(obj)mean(obj.(param)),inRegion(arrayfun(@(obj)length(obj.(param)),inRegion)>=1));
                [values, edges] = histcounts(meanParam,100,'Normalization','Probability','BinLimits',[1,100]);
                h = plot(ax,edges(2:end),values,'-','Color',colorOrder(n,:));
            end
            hold(ax,'off')
        end
        function h = plotDistribution(obj,ax,varargin)
            if isempty(obj.regions)
                obj = obj.setupRegions;
            end
            if ismember('volume',varargin(cellfun(@(c)isa(c,'char'),varargin)))
                param = 'volume';
            else
                param = 'size';
            end
            if any(strcmpi(varargin,'Vavg')) && strcmp(param,'volume')
                adjusted = @(x,at)at.*x/trapz(at,at.*x);
            elseif any(strcmpi(varargin,'Vavg')) && strcmp(param,'size')
                Fadj = @(r)1e-3*(4/3*pi*r.^3 - 2*(pi/3*max([r-obj.scale.depth/2;zeros(size(r))],[],1).^2.*(2*r+obj.scale.depth/2)));
                adjusted = @(x,at)x.*Fadj(at)/trapz(at,x.*Fadj(at));
            else
                adjusted = @(x,at)x;
            end
            
            if isempty(varargin(cellfun(@(c)isa(c,'numeric'),varargin)))
                regionIDs = 1:numel(obj.regions);
            else
                regionIDs = [varargin{cellfun(@(c)isa(c,'numeric'),varargin)}];
            end
            
            if isempty(regionIDs)
                return
            else
                objects = obj.Objects;
            end
            
            if any(strcmpi(varargin,'removeIdle'))
                objects = videoMaker.removeStationary(objects);
            end
            
            if nargin < 2 || isempty(ax)
                ax = axes(figure);
            end
            
            % filter objects based on location of droplets
            hold(ax,'on')
            colorOrder = get(ax, 'ColorOrder');
            for n = regionIDs
                if isempty(obj.rect)
                    regionXs = obj.regions{n}(1:2:end-1); regionYs = obj.regions{n}(2:2:end);
                else
                    regionXs = obj.regions{n}(1:2:end-1) - obj.rect(1); regionYs = obj.regions{n}(2:2:end) - obj.rect(2);
                end
                % inRegion = objects(cellfun(@(pos)any(inpolygon(pos(:,1),pos(:,2),regionXs,regionYs)),{objects.position}));
                inRegion = @(object)inpolygon(object.position(:,1),object.position(:,2),regionXs,regionYs);
                
%                 if n == 1 % added for AG1B-5
%                     inRegion(arrayfun(@(obj)mean(obj.size),inRegion)<50) = [];
%                 end
                meanParam = arrayfun(@(object)mean(object.(param)(inRegion(object))),objects(arrayfun(@(object)sum(inRegion(object))>0,objects)));
                %meanParam = cell2mat(arrayfun(@(object)object.(param)(inRegion(object)),objects(arrayfun(@(object)sum(inRegion(object))>2,objects)),'uni',0));
                %2 meanParam = cell2mat(arrayfun(@(obj)(obj.(param)),inRegion(arrayfun(@(obj)length(obj.(param)),inRegion)>=1),'uni',0));
                % meanParam = arrayfun(@(obj)mean(obj.(param)),inRegion(arrayfun(@(obj)length(obj.(param)),inRegion)>=1));
                [values, at] = ksdensity(meanParam,sort([linspace(0,100,100),meanParam']));
                h = plot(ax,at,adjusted(values,at),'-','Color',colorOrder(n,:));
            end
            hold(ax,'off')
        end
        function h = plotFrame(obj,frameID)
            if frameID > length(obj.frames)
                return
            end
            v = obj.openVideo;
            frameToShow = obj.frames(frameID);
            v.CurrentTime = frameToShow.time;
            
            frameData = frameToShow.reportData;
            centers = frameData(:,1:2) + repmat(obj.rect(1:2),size(frameData,1),1); radii = frameData(:,3);
            
            h = axes(figure); imshow(readFrame(v)); viscircles(centers, radii);
        end
        function plotRegions(obj,ax,varargin) 
            if isempty(obj.regions)
                msgbox('There are currently no regions, use "obj.setupRegions" to add regions.');
                return
            end
                 
            if isempty(varargin(cellfun(@(c)isa(c,'numeric'),varargin)))
                regionIDs = 1:numel(obj.regions);
            else
                regionIDs = varargin{1};
            end
            
            if isempty(regionIDs)
                return
            end
            
            if nargin < 2 || isempty(ax)
                % Get example image
                v = obj.openVideo;
                v.CurrentTime = v.Duration/2;
                Im = readFrame(v);
                ax = axes(figure);
                imshow(Im);
            end
            
            % filter objects based on location of droplets
            hold(ax,'on')
            colorOrder = get(ax, 'ColorOrder');
            for n = regionIDs
                regionXs = obj.regions{n}(1:2:end-1); regionYs = obj.regions{n}(2:2:end);
                
                plot(ax,regionXs,regionYs,'-','LineWidth',2,'Color',colorOrder(n,:));
                text(ax,max(regionXs+7),mean(regionYs),num2str(n),'FontSize',16,'FontWeight','bold','Color',colorOrder(n,:));
            end
            hold(ax,'off')
        end
        
        % Dependent variable functions
        function objects = get.Objects(obj)
            if ~isfield(obj.analysis,'tracking')
                objects = [];
                return
            end
            
            % Reorient data for tracking
            frameData = arrayfun(@(n)reportData(obj.frames(n)),1:numel(obj.frames),'uni',0);
            frameData = cellfun(@(x,y)[x repmat(y,size(x,1),1)],frameData,num2cell(1:numel(frameData)),'uni',0);
            particlePositions = sortrows(cell2mat(reshape(frameData,[],1)),[4 1 2]);
            % particlePositions(:,1:2) = particlePositions(:,1:2) + repmat(obj.rect(1:2),size(particlePositions,1),1);
            
            % Retrieve particle tracking data and sort
            particleTracks = obj.analysis.tracking;
            particleIDs = sortrows(particleTracks,[3 1 2]);
            
            particles = splitapply(@(x){x},particlePositions,particleIDs(:,end));
            vol = @(r)1e-3*(4/3*pi*r.^3 - 2*(pi/3*max([r-obj.scale.depth/2,zeros(size(r))],[],2).^2.*(2*r+obj.scale.depth/2)));
            
            objects = cellfun(@(P)struct('size',P(:,3)*obj.scale.factor,'position',P(:,[1,2]),...
                                         'volume', vol(P(:,3)*obj.scale.factor),...
                                           'atTime',P(:,4),'inFrame',P(:,5)),particles);
        end
    end
    %% Hidden methods
    methods (Hidden = true) 
        % Video access functions
        function v = openVideo(obj)
            try v = VideoReader(obj.videoLink);
            catch me
                msgbox({'Error while opening video link:' me.message})
            end
        end
        function obj = prepareForOdyssee(obj)
            [~, name] = fileparts(obj.videoLink);
            odysseeFolder = '/Users/Soren/Desktop'; %'/Volumes/homes/home00/sbrandt/Matlab';
            if ~exist(odysseeFolder,'dir')
                error('Odyssee folder not found.')
            end
            if exist([odysseeFolder,'/',name],'dir')
                error('Video folder already exists.')
            else
                mkdir([odysseeFolder,'/',name]);
            end
            
            video = obj.openVideo;
            estFrames = ceil(video.Duration*video.FrameRate);

            %---------------------------------------------------------------
            %         Step 1 Initialize frames
            %---------------------------------------------------------------
            if isempty(obj.imageType)
                obj.imageType = 'fluorescence'; %videoMaker.determineImageType(firstImage);
            end
            
            %---------------------------------------------------------------
            %         Step 2 Continue for all frames
            %---------------------------------------------------------------
            progress = waitbar(0,'initializing frames');

            Nr = 1; timeRem = nan(1,100);
            while hasFrame(video) && video.CurrentTime < video.Duration*1.05
                tic % start clock
                currentTime = video.CurrentTime;
                waitbar(currentTime/video.Duration,progress,...
                    ['analyzing frame ',num2str(Nr),' of ',num2str(estFrames),...
                    ' (time rem: ',num2str(round((estFrames-Nr)*nanmean(timeRem)/60)),'min)']);
                             
                currentImage = readFrame(video);
                
                imwrite(currentImage,[odysseeFolder,'/',name,'/',num2str(currentTime),'.tif'])
                
                Nr = Nr +1;
                timeRem(rem(Nr,100)+1) = toc;
            end
            delete(progress)
            
            % save copy to odyssee
            obj.save([odysseeFolder,'/',name]);
        end
        
        % Video analysis functions
        function obj = setupVideoAnalysisFor(obj,varargin)
            % Setup frame analysis
            implementedTypes = {'circle','full'};
            obj.forType = implementedTypes(ismember(implementedTypes,lower(varargin)));
            
            implementedImageTypes = {'fluorescence', 'brightfield', 'darkfield'};
            obj.imageType = [implementedImageTypes{ismember(implementedImageTypes,lower(varargin))}];
            
            % Setup analysis
            implementedAnalysis = {'tracking','distribution','streamlines'};
            obj.forAnalysis = implementedAnalysis(ismember(implementedAnalysis,lower(varargin)));
            
            % Setup rectangle
            if any(ismember(varargin,{'rect', 'rescale'}))
                v = obj.openVideo;
                v.CurrentTime = v.Duration/2;
                
                h = figure;
                [~,obj.rect] = imcrop(readFrame(v));
                delete(h)
            end
            
            % Setup regions
            if any(strcmpi(varargin,'regions'))
                obj = obj.setupRegions;
            end
            
            % Setup scale
            if isempty(obj.scale)
                obj = obj.setupScale;
            end
            
            % Return remainder of varargin as options
            lookFor = [implementedTypes, implementedImageTypes, implementedAnalysis, {'rect', 'regions'}];
            obj.options = varargin(~ismember(varargin,lookFor));
        end
        function frames = analyseFrames(obj)
            video = obj.openVideo;
            estFrames = ceil(video.Duration*video.FrameRate);
            if any(contains(obj.options,{'show:'},'IgnoreCase',true))
                showFreq = str2double(cell2mat(...
                    regexp(obj.options{contains(obj.options,{'show:'},'IgnoreCase',true)},'\d*','Match')));
            else
                showFreq = inf;
            end
            
            %---------------------------------------------------------------
            %         Step 1 Initialize frames
            %---------------------------------------------------------------
            progress = waitbar(0,'initializing frames');
            firstImage = readFrame(video);
            clearFeatures = false(size(firstImage)); maxBrightness = max(max(rgb2gray(firstImage)));
            
            if isempty(obj.imageType)
                obj.imageType = videoMaker.determineImageType(firstImage);
            end
            if strcmp(obj.imageType,'brightfield')
                % Pick 5 images throughout the video and multiply binarized
                % image complement
                clearFeatures = true(size(firstImage,1),size(firstImage,2));
                for n = 1:5
                    video.currentTime = video.Duration*n/6;
                    clearFeatures = clearFeatures.*imbinarize(imcomplement(rgb2gray(readFrame(video))));
                end
                clearFeatures = repmat(bwareafilt(logical(clearFeatures),1),1,1,size(firstImage,3));
                
                % Reset video to second frame
                video.currentTime = 0; readFrame(video);
                firstImage(clearFeatures) = maxBrightness;
            end
            
            currentFrame = videoFrame(firstImage, obj.rect, 0, obj.forType{:}, obj.imageType, obj.options{:});
            frames(estFrames) = clean(currentFrame);
            frames = fliplr(frames);
            
            %---------------------------------------------------------------
            %         Step 2 Continue for all frames
            %---------------------------------------------------------------
            Nr = 2; timeRem = nan(1,100);
            while hasFrame(video) && video.CurrentTime < video.Duration*1.05
                tic % start clock
                currentTime = video.CurrentTime;
                waitbar(currentTime/video.Duration,progress,...
                    ['analyzing frame ',num2str(Nr),' of ',num2str(estFrames),...
                    ' (time rem: ',num2str(round((estFrames-Nr)*nanmean(timeRem)/60)),'min)']);
                
                currentImage = readFrame(video);
                currentImage(clearFeatures) = maxBrightness;
                
                currentFrame = videoFrame(currentImage, obj.rect, currentTime, obj.forType{:}, obj.imageType, obj.options{:});
                frames(Nr) = clean(currentFrame);
                
                if rem(Nr,showFreq) == 0
                  plot(currentFrame);
                end
                
                Nr = Nr +1;
                timeRem(rem(Nr,100)+1) = toc;
            end
            
            % clean out empty frames
            frames(arrayfun(@(frame)isempty(frame.analysis),frames)) = [];
            
            delete(progress)
        end
        function obj = setupRegions(obj)
            warning('off'); 
            % Setup regions
            v = obj.openVideo;
            v.CurrentTime = v.Duration/2;
            Im = readFrame(v);

            obj.regions = {};
            done = false;
            while ~done
                h = figure; imshow(Im); title('Close figure to end entering regions')
                [~,posX,posY] = roipoly;

                done = isempty(posX);
                if done
                    break
                else
                    if ~isempty(obj.rect)
                    % adjust ROI to rect
                        R = bbox2points(obj.rect);
                        [posX, posY] = polybool('intersection', posX, posY, R(:,1), R(:,2));
                    end

                    % store and show
                    close(h)
                    obj.regions{end+1} = reshape([posX';posY'],1,[]);
                    Im = insertShape(Im,'Polygon',obj.regions);
                end
            end
            warning('on');
        end
        function obj = setupScale(obj)
            % Setup scale
            v = obj.openVideo;
            v.CurrentTime = v.Duration/2;
            Im = readFrame(v);
            
            fig = figure; imshow(Im); title('Select scale')
            
            try [x,y,~] = improfile;
            catch
                datacursormode toggle
                return
            end
            close(fig)
            distance = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
            
            input = inputdlg({['Channel width /mm (',num2str(distance),'pixels)'],'Channel depth /um'},'Setup scale',1,{'1.5','25'});
            obj.scale = struct('factor',str2double(input{1})*1000/distance,'depth',str2double(input{2}));
        end
        
        % Data analysis functions
        function analysis = analyseVideo(obj)
            %---------------------------------------------------------------
            %         Step 0 Determine analysis to be performed
            %---------------------------------------------------------------
            %ie if streamlines, also perform tracking
            
            %---------------------------------------------------------------
            %         Step 1 Perform analysis
            %---------------------------------------------------------------
            analysis = obj.analysis;
            for n = 1:numel(obj.forAnalysis)
                of = obj.forAnalysis{n};
                analysis.(of) = runAnalysis(obj,of);
            end
        end
        function analysis = trackParticles(obj)
            % analysis = trackParticles(obj) performs particle tracking
            % analysis on the objects found by videoMaker. Use runAnalysis
            % or redoAnalysis to perform trackparticles with options.
            %
            %   analysis = trackParticles(obj) returns tracked, indexed
            %   positions. See track for details.
            %
            %   obj = redoAnalysis(obj,'tracking') performs particle
            %   tracking analaysis for obj.
            %
            %   obj = redoAnalysis(obj,'tracking','mem: 2','maxdisp: 5') 
            %   performs particle tracking with options. 
            %       'mem: n' adjusts particle memory, where particle 
            %       positions can disappear for n time steps before the 
            %       particle is considered gone.
            %       'maxdisp: d' adjusts the maximum displacement of the
            %       particles. d should be on the order of the particles as
            %       a first estimate. 
            %       See track for more details.
            %
            %   See also videoMaker.runAnalysis, videoMaker.trackParticles.
            
            % Reorient data for tracking
            frameData = arrayfun(@(n)reportData(obj.frames(n)),1:numel(obj.frames),'uni',0);
            particlePositions = cell2mat(reshape(frameData,[],1));
            
            % set memory and displacement
            if any(contains(obj.options,{'mem:'},'IgnoreCase',true))
                mem = str2double(cell2mat(...
                    regexp(obj.options{contains(obj.options,{'mem:'},'IgnoreCase',true)},'\d*','Match')));
            else; mem = 2; 
            end
            if any(contains(obj.options,{'maxdisp:'},'IgnoreCase',true))
                maxdisp = str2double(cell2mat(...
                    regexp(obj.options{contains(obj.options,{'maxdisp:'},'IgnoreCase',true)},'\d*','Match')));
            else; maxdisp = 5; 
            end
            
            % track particle positions
            % addPath('tracking');
            param = struct('mem',mem,'good',0,'dim',2,'quiet',0);
            analysis = track(particlePositions(:,[1 2 4]),maxdisp,param); % median(particlePositions(:,3))
            %uniqueTracks = unique(res(:,4));
            %tracks = arrayfun(@(n)struct('x',res(res(:,4) == n,1),'y',res(res(:,4) == n,2),'time',res(res(:,4) == n,3)),uniqueTracks);
            
            % return analysis
            %analysis = struct('for','tracking',...
                        %'tracks',tracks);
        end
        function obj = sizeDistributions(obj) % WRITE include regions
        end
    end
    %% Static, hidden functions
    methods (Static, Hidden = true)
        % Static or undefined functions
        function plotTracks(obj,ax,frame) % Remove access to obj and inlcude in video output
            trail = 0.5; % second
            % check if analysis contains else return warning and do nothing
            if ~isfield(obj.analysis,'tracking')
                warning('No tracking analysis found.'); return
            end
            
            % make sure arguments are populated
            if nargin < 2 || isempty(ax)
                ax = gca;
            end
            if nargin < 3 || frame > length(obj.frames)
                time = inf;
            else
                time = obj.frames(frame).time;
            end
            
            % Plot tracking results into axis handle
            trackingResults = obj.analysis.tracking;
            current = trackingResults(time > trackingResults(:,3) & trackingResults(:,3) > time-trail,:);
            [~, ~, trackMap] = unique(current(:,4));
            tracks = splitapply(@(x){x},current(:,[1,2]),trackMap);
            XYs = cellfun(@(track){track(:,1), track(:,2)},tracks,'uni',0); XYs = [XYs{:}];
            
%             current = cellfun(@(x)any(x < time & x > time-trail),{obj.analysis.tracking.tracks.time});
%             times = {obj.analysis.tracking.tracks(current).time};
%             tracks = [{obj.analysis.tracking.tracks(current).x};{obj.analysis.tracking.tracks(current).y}];
%             for n = 1:sum(current)
%                 remove = times{n} >= time | times{n} <= time-trail;
%                 tracks{1,n}(remove) = []; tracks{2,n}(remove) = [];
%             end
            if ~isempty(tracks)
                hold(ax,'on')
                plot(ax,XYs{:},'Color','r')
                hold(ax,'off')
            end
        end
        function Im = addTracksToFrame(obj,frame,Im,color) % Remove access to obj and inlcude in video output
            trail = 10; % second
            % check if analysis contains else return warning and do nothing
            if ~isfield(obj.analysis,'tracking')
                warning('No tracking analysis found.'); return
            end
            
            % make sure arguments are populated
            if nargin < 2 || frame > length(obj.frames)
                time = inf;
            else
                time = obj.frames(frame).time;
            end
            if nargin < 4
                color = [1 1 1]*255;
            else
                color = color*255;
            end
            
            % Plot tracking results into axis handle
            trackingResults = obj.analysis.tracking;
            current = trackingResults(time > trackingResults(:,3) & trackingResults(:,3) > time-trail,:);
            trackPoints = splitapply(@(x){round(x)},current(:,[1,2]),current(:,4));
            
%             current = cellfun(@(x)any(x < time & x > time-trail),{obj.analysis.tracking.tracks.time});
%             times = {obj.analysis.tracking.tracks(current).time};
%             tracks = {obj.analysis.tracking.tracks(current).x; obj.analysis.tracking.tracks(current).y};
%             for n = 1:sum(current)
%                 remove = times{n} >= time | times{n} <= time-trail;
%                 tracks{1,n}(remove) = []; tracks{2,n}(remove) = [];
%             end
            
            % trackPoints = arrayfun(@(n)unique(round(horzcat(tracks{1:2,n})),'rows'),1:size(tracks,2),'uni',0);
            trackDiffs = cellfun(@(x)max(abs(diff(x,1,1))+1,[],2),trackPoints,'uni',0);
            for m = find(cellfun(@(x)any(x > 2),trackDiffs)) %1:size(trackDiffs,2)
                addPoints = trackDiffs{m} > 2; nrOfPoints = trackDiffs{m}(addPoints);
                newPointsX = cell2mat(arrayfun(@(a,b,c)round(linspace(a,b,c)),trackPoints{m}([addPoints;false],1),trackPoints{m}([false;addPoints],1),nrOfPoints,'uni',0)')';
                newPointsY = cell2mat(arrayfun(@(a,b,c)round(linspace(a,b,c)),trackPoints{m}([addPoints;false],2),trackPoints{m}([false;addPoints],2),nrOfPoints,'uni',0)')';
                trackPoints{m} = [trackPoints{m};[newPointsX newPointsY]];
            end
            markPoints = cell2mat(reshape(trackPoints,[],1));
            if ~isempty(markPoints)
                mark = false(size(Im(:,:,1)));
                mark(sub2ind(size(mark),markPoints(:,2),markPoints(:,1))) = 1;
                A = Im(:,:,1); A(mark) = color(1);
                B = Im(:,:,2); B(mark) = color(2);
                C = Im(:,:,3); C(mark) = color(3);
                Im(:,:,1) = A;
                Im(:,:,2) = B;
                Im(:,:,3) = C;
            end
        end
        function type = determineImageType(Im)
            if nargin < 1
                type = []; return
            elseif size(Im,3) > 1
                Im = rgb2gray(Im);
            end
            
            binaryIm = imbinarize(Im);
            ImROI = bwpropfilt(binaryIm,'Eccentricity',[0 0.5]);
            binaryImROI = bwconvhull(ImROI); %regionprops(ImROI,'BoundingBox'); binaryImROI = binaryImROI.BoundingBox;
            
            inROI = mean(mean(binaryIm(ImROI)));
            outsideROI = mean(mean(binaryIm(~ImROI & binaryImROI)));
            
%             binaryImComplement = imcomplement(binaryIm);
%             binaryImBackground = imfill(edge(binaryIm),'holes');
%             binaryImComplementBackground = imfill(edge(binaryImComplement),'holes');
%             
%             ImVariation = std(binaryImBackground(:));
%             ImComplementVariation = std(binaryImComplementBackground(:));
            
            if round(abs(inROI - outsideROI),1) < 0.5
                type = 'brightfield';
            elseif inROI < outsideROI
                type = 'darkfield'; % not sure about this one
            else
                type = 'fluorescence';
            end
        end
        function objects = removeStationary(objects)
            if length(objects) < 1
                return
            end
            
            % find objects that move out of defined radius during the video
            distance = @(O,P)sqrt( (P(:,1) - O(1)).^2 + (P(:,2) - O(2)).^2 );
            
            kept = arrayfun(@(object)any( distance(object.position(1,:),object.position) > object.size(1) ),objects);
            
            % keep valuable objects (i.e. moving)
            objects = objects(kept);
        end
    end
end
