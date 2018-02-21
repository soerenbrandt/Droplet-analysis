classdef videoMaker
    % Analysis of objects in video and returns data or video with
    % respective features and analysis.
    %
    % Object detection (based on videoFrame class) includes
    %   'circles'                     - detects circular objects in video frame
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
    end
    properties (Dependent)
        Distribution % distribution of particles in video
    end
    properties (SetAccess = protected)
        forType % analysis performed for object types
        forAnalysis % analysis methods performed
        rect % rectangle to crop the video
        regions % relevant regions in image
    end
    properties (SetAccess = immutable)
        date % creation date of the video file
    end
    
    methods
        function obj = videoMaker(videoLink,varargin) %WRITE
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
            
            %---------------------------------------------------------------
            %         Step 2 Set analysis variables
            %---------------------------------------------------------------
            obj = setupVideoAnalysisFor(obj,varargin{:});
            
            %---------------------------------------------------------------
            %         Step 3 Analyse frames
            %---------------------------------------------------------------
            obj.frames = analyseFrames(obj);
            
            %---------------------------------------------------------------
            %         Step 4 Analyse video
            %---------------------------------------------------------------
            obj.analysis = analyseVideo(obj);
            
        end
        function obj = redoAnalysis(obj,varargin)
            %---------------------------------------------------------------
            %         Step 1 Set analysis variables
            %---------------------------------------------------------------
            obj = setupVideoAnalysisFor(obj,varargin);
            
            %---------------------------------------------------------------
            %         Step 2 Analyse frames
            %---------------------------------------------------------------
            obj.frames = analyseFrames(obj);
            
            %---------------------------------------------------------------
            %         Step 3 Analyse video
            %---------------------------------------------------------------
            obj.analysis = analyseVideo(obj);
        end
        function save(obj)
            [path,name] = fileparts(obj.videoLink);
            
            save([path,'/',name,'.mat'],'obj')
        end
        function analysis = runAnalysis(obj,of)
            switch of
                case 'tracking'
                    analysis = trackParticles(obj);
                otherwise
                    analysis = [];
                    return
            end
        end
        function makeVideo(obj,varargin)
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
            video = videoMaker.openVideo(obj);
            
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
        function makeCVS(obj,varargin) %WRITE
        end
	end
    methods (Hidden = true)
        function frames = analyseFrames(obj)
            video = videoMaker.openVideo(obj);
            estFrames = ceil(video.Duration*video.FrameRate);
            
            progress = waitbar(0,'initializing frames');
            currentFrame = videoFrame(readFrame(video), obj.rect, 0, obj.forType{:});
            frames(estFrames) = clean(currentFrame);
            frames = fliplr(frames);
            Nr = 2; timeRem = nan(1,estFrames);
            while hasFrame(video) && video.CurrentTime < video.Duration*1.05
                tic % start clock
                currentTime = video.CurrentTime;
                waitbar(currentTime/video.Duration,progress,...
                    ['analyzing frame ',num2str(Nr),' of ',num2str(estFrames),...
                    ' (time rem: ',num2str(round((estFrames-Nr)*nanmean(timeRem)/60)),'min)']);
                
                currentFrame = videoFrame(readFrame(video), obj.rect, currentTime, obj.forType{:});
                frames(Nr) = clean(currentFrame);
                
                Nr = Nr +1;
                timeRem(Nr) = toc;
            end
            
            delete(progress)
        end
        function obj = setupVideoAnalysisFor(obj,varargin)
            % Setup frame analysis
            implementedTypes = {'circle'};
            obj.forType = implementedTypes(ismember(implementedTypes,lower(varargin)));
            
            % Setup analysis
            implementedAnalysis = {'tracking','distribution','streamlines'};
            obj.forAnalysis = implementedAnalysis(ismember(implementedAnalysis,lower(varargin)));
            
            % Setup rectangle
            if any(strcmpi(varargin,'rect'))
                v = videoMaker.openVideo(obj);
                v.CurrentTime = v.Duration/2;
                
                h = figure;
                [~,obj.rect] = imcrop(readFrame(v));
                delete(h)
            end
            
            % Setup regions
            if any(strcmpi(varargin,'regions'))
                v = videoMaker.openVideo(obj);
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
                        close(h)
                        obj.regions{end+1} = reshape([posX';posY'],1,[]);
                        Im = insertShape(Im,'Polygon',obj.regions);
                    end
                end
            end
        end
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
                of = obj.forAnalysis(n);
                analysis.(of) = runAnalysis(obj,of);
            end
        end
        function analysis = trackParticles(obj)
            % Reorient data for tracking
            frameData = arrayfun(@(n)report(obj.frames(n)),1:numel(obj.frames),'uni',0);
            particlePositions = cell2mat(reshape(frameData,[],1));
            
            % track particle positions
            res = track(particlePositions,10);
            uniqueTracks = unique(res(:,4));
            tracks = arrayfun(@(n)struct('x',res(res(:,4) == n,1),'y',res(res(:,4) == n,2),'time',res(res(:,4) == n,3)),uniqueTracks);
            
            % return analysis
            analysis = struct('for','tracking',...
                        'tracks',tracks);
        end % WRITE
        function obj = sizeDistributions(obj) % WRITE include regions
        end
    end
    methods (Static, Hidden = true)
        function v = openVideo(obj)
            try v = VideoReader(obj.videoLink);
            catch me
                msgbox({'Error while opening video link:' me.message})
            end
        end
        function plotTracks(obj,ax,frame)
            trail = 10; % second
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
            current = cellfun(@(x)any(x < time & x > time-trail),{obj.analysis.tracking.tracks.time});
            times = {obj.analysis.tracking.tracks(current).time};
            tracks = [{obj.analysis.tracking.tracks(current).x};{obj.analysis.tracking.tracks(current).y}];
            for n = 1:sum(current)
                remove = times{n} >= time | times{n} <= time-trail;
                tracks{1,n}(remove) = []; tracks{2,n}(remove) = [];
            end
            if ~isempty(tracks)
                hold(ax,'on')
                plot(ax,tracks{:},'Color','w')
                hold(ax,'off')
            end
        end
        function Im = addTracksToFrame(obj,frame,Im,color)
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
            current = cellfun(@(x)any(x < time & x > time-trail),{obj.analysis.tracking.tracks.time});
            times = {obj.analysis.tracking.tracks(current).time};
            tracks = {obj.analysis.tracking.tracks(current).x; obj.analysis.tracking.tracks(current).y};
            for n = 1:sum(current)
                remove = times{n} >= time | times{n} <= time-trail;
                tracks{1,n}(remove) = []; tracks{2,n}(remove) = [];
            end
            trackPoints = arrayfun(@(n)unique(round(horzcat(tracks{1:2,n})),'rows'),1:size(tracks,2),'uni',0);
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
        function plotDistribution(obj,ax) %WRITE
        end
        function plotRegions(obj,ax) %WRITE
        end
    end
end
