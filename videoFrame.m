classdef videoFrame
    % Comes with VittoPowerAnalysisV2 to interpret droplet shape and size
    properties
        % might not be worth it... stats % regionprops of the image after image processing
        analysis % results based on shape request
        time % currentTime in video of frame
    end
    properties (Dependent)
        Positions % Locations of objects in image
        Radii % radius dimension of the objects (in px)
    end
    properties (SetAccess = protected)
        type % Image type (fluorescence or brightfield)
        rect % Image rectangle
        Image % original image given
    end
    
    methods
        function obj = videoFrame(Im, rect, time, varargin) % finds approximate region properties in image for further analysis
            %---------------------------------------------------------------
            %         Step 1 Initialize immutable properties
            %---------------------------------------------------------------
            if nargin == 0
                return
            end
            
            if nargin > 1 && ~isempty(rect)
                obj.rect = rect;
                Im = imcrop(Im,rect);
            end
            if nargin > 2
                obj.time = time;
            end
            obj.Image = rgb2gray(Im);
            
            if any(cellfun(@(str)strcmp(str,'brightfield'),varargin))
                obj.type = 'brightfield';
            else
                obj.type = 'fluorescence';
            end
            
            %---------------------------------------------------------------
            %         Step 2 Aproximate properties of image
            %           (skipping this for now since it is slow)
            %---------------------------------------------------------------
            % obj.stats = approximateFrameProperties(obj,'Area','WeightedCentroid','Eccentricity','MajorAxisLength','MinorAxisLength');
            
            %---------------------------------------------------------------
            %         Step 3 Run analysis
            %---------------------------------------------------------------
            implementedAnalysis = {'circle'}; % potentially ellipse, droplet, rectangle
            performAnalysisFor = find(ismember(implementedAnalysis,lower(varargin)));
            
            for n = performAnalysisFor
                forType = implementedAnalysis{n};
                obj.analysis.(forType) = runAnalysis(obj,forType);
            end
        end
        function obj = clean(obj)
            obj.Image = [];
        end
        function h = plot(obj,ax,varargin)
%             if isempty(obj.analysis)
%                 h = plot(obj.stats.WeightedCentroid(1), obj.stats.WeightedCentroid(2),'ro');
%                 return
%             end
            
            if nargin < 2 || isempty(ax)
                ax = axes(figure);
                imshow(obj.Image)
            end
            
            % CHANGE if varargin specifies anything, else plot first item
            if isempty(obj.analysis)
                return
            else
                hold(ax,'on')
                analysisFor = fieldnames(obj.analysis);
                switch analysisFor{1}
                    case 'circle'
                        if any(cellfun(@(str)strcmp(str,'noOutline'),varargin))
                            h = plot(obj.analysis.circle.centers(:,1),obj.analysis.circle.centers(:,2),'r+');
                        else
                            h = viscircles(ax,obj.analysis.circle.centers,obj.analysis.circle.radii);
                        end
                    otherwise
                        return
                end
            end
            hold(ax,'off')
        end
        function analysis = runAnalysis(obj,type)
            switch type
                case 'circle'
                    if strcmp(obj.type,'brightfield')
                        analysis = detectBrightfieldDroplets(obj);
                    else
                        analysis = detectDroplets(obj);
                    end
                otherwise
                    analysis = [];
                    return
            end
        end
        function report = report(obj,type)
            analysisFor = fieldnames(obj.analysis);
            if nargin < 2
                type = analysisFor{1};
            elseif ~ismember(analysisFor, type)
                report = []; return
            end
            
            report = obj.analysis.(type).centers;
            report(:,3) = obj.time;
        end
	end
    methods (Hidden = true)
        function analysis = detectDroplets(obj)
            % setJump = 3;
            warning('off')
            
            % TRY AGAIN ISOLATING DROPLETS AND FITTING ONLY RADII THAT FIT
            % INTO SIZE OF REGION
            %---------------------------------------------------------------
            %         Step 1 Isolate droplet regions
            %---------------------------------------------------------------
%             A = obj.Image > 30; % thresholding original image
%             B = bwareaopen(A,15); % excluding small regions
%             C = regionprops(A,'BoundingBox'); % isolating droplet regions
%             
%             image = cell(uint8(zeros(size(obj.Image))));
%             for n = 1:numel(C)
%                 image{n} = imcrop(obj.Image,C(n).BoundingBox);
%             end
            
            %images{1} = obj.Image; images{2} = imresize(obj.Image,2);
            rescale = table([0.93; 0.93; 0.93],[10 30;4 10;1 4],'VariableNames',{'Scale','inRange'});
            
            %---------------------------------------------------------------
            %         Step 1 Set radius relations (Deprecated: Approximate radii from stats)
            %---------------------------------------------------------------
%             circularObjects = [obj.stats.Eccentricity] <= 0.75;
%             radiusRange = [min([1 floor([obj.stats(circularObjects).MajorAxisLength]/4 + [obj.stats(circularObjects).MinorAxisLength]/4)]) ...
%                            max(ceil([obj.stats(circularObjects).MajorAxisLength]/4 + [obj.stats(circularObjects).MinorAxisLength]/4))];
            
            % range = radiusRange(1)*[1 setJump];
            %---------------------------------------------------------------
            %         Step 2 Find circles in ranges
            %---------------------------------------------------------------
            centers = nan(100,2); radii = nan(100,1); metric = nan(100,1);
            
            % while max(range) < radiusRange(2)
            for n = 1:height(rescale)
                sensitivity = rescale{n,1}; inRange = rescale{n,2};
                [C, R, M] = videoFrame.detectCircles(obj.Image,inRange,sensitivity);%0.92);
                
                % rescale center and radius
                %C = C./scale; R = R./scale;
                
                currentIndex = max([find(radii(~isnan(radii)),1,'last'),0])+1;
                
                for m = 1:numel(R)
                    overlap = (C(m,1) - centers(:,1)).^2 + (C(m,2) - centers(:,2)).^2 < radii.^2;
                    
                    if any(overlap)
                        % replace better fits
                        centers(overlap & metric < M(m),:) = repmat(C(m,:),sum(overlap & metric < M(m)),1);
                        radii(overlap & metric < M(m),:) = R(m);
                        metric(overlap & metric < M(m),:) = M(m);
                    else
                        % add new circle
                        centers(currentIndex,:) = C(m,:);
                        radii(currentIndex) = R(m);
                        metric(currentIndex) = M(m);
                        currentIndex = currentIndex + 1;
                    end
                end
                
                % range = range(2)*[1 setJump];
            end
            
            warning('on')
            % Clean up centers and radii (remove identical doubles)
            results = unique([centers(metric > 0.1,:), radii(metric > 0.1)],'rows');
            analysis = struct('for','circles',...
                        'centers',results(~isnan(results(:,3)),1:2),...
                        'radii',results(~isnan(results(:,3)),3));
        end
        function analysis = detectBrightfieldDroplets(obj)
            A = imbinarize(imcomplement(obj.Image)); % binarize the inverse brightfield image
            B = imfill(A,'holes'); % fill any fully enclosed white regions
            C = B; C(A) = 0; % remove outlines of said regions
            stats = regionprops(C, obj.Image, 'Area','WeightedCentroid','Eccentricity','MajorAxisLength','MinorAxisLength');
            stats([stats.Eccentricity] >= 0.75) = [];
            D = xor(bwareaopen(A,1), bwareaopen(A,10*max([stats.Area]))); % isolate dark droplet rim in original image
            d = bwdist(imcomplement(D)); d = median(d(d>0)); % d is the half-width of the dark rim of the droplet
            
            analysis = struct('for','circles',...
                        'centers',vertcat(stats.WeightedCentroid),...
                        'radii',mean([vertcat(stats.MajorAxisLength),vertcat(stats.MinorAxisLength)]/2,2) +2*d);
        end
        function stats = approximateFrameProperties(obj,varargin) % WORK IN PROGRESS (TOO SLOW)
            A = obj.Image;
            A(A < 30) = 0; % thresholds low intensity parts of image
%             B = imbinarize(A,0.6);
%             C = imcomplement(B);
%             tic; D = C; for n = 1:9; D = imclose(D,strel('line',7,20*n)); end; toc
%             F = bwmorph(B,'hbreak',inf);
            B = imtophat(imgaussfilt(A),strel('disk',15)); % smooths image and removes non-uniform brightness
            C = imopen(localcontrast(B,1,1),strel('disk',1)); % increases contrast locally and smooth single pixel objects
            D = imbinarize(C,0.6); % creates black-and-white image
            E = bwmorph(D,'hbreak',inf); % removes small connections between regions
            F = imfill(E,'holes'); % fills in holes in regions (from local contrast algorithm)
            stats = regionprops(F,A,varargin{:});
        end
    end
    methods (Static, Hidden = true)
        function [centers, radii, metric] = detectCircles(Im,radius,varargin)
            % finds circles in Im

            % validate input
            validateattributes(Im, {'numeric'},{'3d'})
            validateattributes(radius, {'numeric'},{'nondecreasing','2d'})
            if length(radius) == 1
                radius = [radius*0.9 radius*1.1];
            else
                validateattributes(radius, {'numeric'},{'size',[1,2]})
            end

            % check if Image is RGB and needs to be converted still
            if size(Im,3) == 3
                Im = rgb2gray(Im);
            end

            % check input parameters
            if strcmp('dark',varargin); polarity = 'dark'; else; polarity = 'bright'; end
            if strcmp('TwoStage',varargin); method = 'TwoStage'; else; method = 'PhaseCode'; end
            sensitivity = [varargin{cellfun(@(y)isa(y,'double'),varargin)}, 0.85];

            % find circles
            [centers, radii, metric] = imfindcircles(Im,radius,'ObjectPolarity',polarity,'Sensitivity',sensitivity(1), 'Method', method);
        end
    end
end