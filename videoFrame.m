classdef videoFrame
    % videoFrame V1.4
    % minor bug fixes
    
    % Comes with videoMaker to analyze droplets moving in image sequence
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
        rescaled % true if image was rescaled for analysis (by 2x)
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
                if any(cellfun(@(str)strcmp(str,'rescale'),varargin))
                    obj.rescaled = true;
                    Im = imcrop(imresize(Im,2),[2*rect(1) 2*rect(2)+rect(4)/2 rect(3)*2 rect(4)]);
                else
                    Im = imcrop(Im,rect);
                end
            end
            if nargin > 2
                obj.time = time;
            end
            if size(Im,3) == 3
                obj.Image = rgb2gray(Im);
            else
                obj.Image = Im;
            end
            
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
            implementedAnalysis = {'circle','full'}; % potentially ellipse, droplet, rectangle
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
                    case 'full'
                        if any(cellfun(@(str)strcmp(str,'noOutline'),varargin))
                            h = plot(obj.analysis.full.centers(:,1),obj.analysis.full.centers(:,2),'r+');
                        else
                            h = viscircles(ax,obj.analysis.full.centers,obj.analysis.full.radii);
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
                        analysis = detectDroplets(obj); % previously detectDroplets(obj)
                    end
                case 'full'
                    analysis = detectDropletsFULL(obj);
                otherwise
                    analysis = [];
                    return
            end
        end
        function report = reportData(obj,type)
            analysisFor = fieldnames(obj.analysis); 
            if nargin < 2
                type = analysisFor{1};
            elseif ~ismember(analysisFor, type)
                report = []; return
            end
            
            if obj.rescaled
                report = [(obj.analysis.(type).centers(:,1) + obj.rect(1)*0)/2, (obj.analysis.(type).centers(:,2) + obj.rect(4)/2)/2, obj.analysis.(type).radii/2];
            else
                report = [obj.analysis.(type).centers, obj.analysis.(type).radii];
            end
            report(:,4) = obj.time;
        end
        function positions = positions(obj,type)
            analysisFor = fieldnames(obj.analysis);
            if nargin < 2
                type = analysisFor{1};
            elseif ~ismember(analysisFor, type)
                positions = []; return
            end
            
            if obj.rescaled
                positions = [(obj.analysis.(type).centers(:,1) + obj.rect(1)*0)/2, (obj.analysis.(type).centers(:,2) + obj.rect(4)/2)/2];
            else
                positions = obj.analysis.(type).centers;
            end
            positions(:,3) = obj.time;
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
            rescale = table([0.93; 0.93; 0.93; 0.93],[30 60;10 30;4 10;1 4],'VariableNames',{'Scale','inRange'});
            if obj.rescaled
                rescale = table(rescale{:,1},rescale{:,2}*2,'VariableNames',{'Scale','inRange'});
            end
            
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
                [C, R, M] = imfindcircles(obj.Image,inRange,'Sensitivity',sensitivity);%0.92);
                
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
        function analysis = detectFluorescentDroplets(obj)
            % threshold image
            A = obj.Image;
            B = adaptthresh(A,0.3); B(B>0.75) = 0.75; % B = A; B(A<25) = 0; B(A>150) = 255; % potentially do without
            C = imbinarize(A,B); % imextendedmax(B,5);

            % find coarse estimate of droplets and filter bad shapes
            % (conjoined droplets)
            imageStats = regionprops(C,A,'BoundingBox','Extent','Image');
            badEstimates = [imageStats.Extent] < 0.75; blank = zeros(size(A)); % start blank guess for droplet locations
            % Note: perfect circles should have extend around 0.78, larger
            % values usually occur for very small droplets

            % filter out conjoined droplets and find a better estimate
            % based on distance maxima
            localRadii = bwdist(~C);
            localMaxima = logical(imregionalmax(localRadii).*bwmorph(imregionalmax(medfilt2(localRadii)),'thicken',2)); % finds maxima and roughly filters local maxima between two droplets
            for i = find(~badEstimates) % remove regions of good estimates from image
                n = imageStats(i);
                Ys = ceil(n.BoundingBox(1)):floor(n.BoundingBox(1)+n.BoundingBox(3));
                Xs = ceil(n.BoundingBox(2)):floor(n.BoundingBox(2)+n.BoundingBox(4));
                localMaxima(Xs,Ys) = zeros(size(n.Image));
                blank(Xs,Ys) = n.Image; % fill up blank guess with good estimates for droplets
            end

            % find located maxima and draw into blank 
            localStats = regionprops(localMaxima,'Centroid','PixelIdxList');
            centers = vertcat(localStats.Centroid);
            radii = arrayfun(@(x)mean(localRadii(x.PixelIdxList))-3,localStats);

            [colums, rows] = meshgrid(1:size(A,2),1:size(A,1));
            inCirc = @(c,r)(colums-c(1)).^2 + (rows-c(2)).^2 < r.^2;

            for n = 1:length(radii)
                blank(inCirc(centers(n,:),radii(n))) = 1;
            end

            % extend the size of droplets to overlap with original image
            finalGuess = bwmorph(blank,'thicken',25);

            F = double(A).*finalGuess; F = F > 255/exp(1); % overlap and threshold guess with original image
            finalStats = regionprops(F, A,'WeightedCentroid','Eccentricity','EquivDiameter');
            finalStats([finalStats.Eccentricity] >= 0.75) = [];

            analysis = struct('for','circles',...
                                    'centers',vertcat(finalStats.WeightedCentroid),...
                                    'radii',vertcat(finalStats.EquivDiameter)/2);
        end
        function analysis = detectBrightfieldDroplets(obj)
            % thresh = double(max(obj.Image(:)) - min(obj.Image(:)))/255 +0.1;
            A = imbinarize(imcomplement(obj.Image)); %,thresh); %0.7); % binarize the inverse brightfield image
            B = imfill(A,'holes'); % fill any fully enclosed white regions
            C = B; C(A) = 0; % remove outlines of said regions
            stats = regionprops(C, obj.Image, 'Area','WeightedCentroid','Eccentricity','MajorAxisLength','MinorAxisLength');
            stats([stats.Eccentricity] >= 0.75) = [];
            D = xor(bwareaopen(A,1), bwareaopen(A,min([10*max([stats.Area]),300000]))); % isolate dark droplet rim in original image
            d = bwdist(imcomplement(D)); d = median(d(d>0)); % d is the half-width of the dark rim of the droplet
            
            analysis = struct('for','circles',...
                        'centers',vertcat(stats.WeightedCentroid),...
                        'radii',mean([vertcat(stats.MajorAxisLength),vertcat(stats.MinorAxisLength)]/2,2) +2*d);
        end
        function analysis = detectDropletsFULL(obj)
            warning('off')
            
            % imfindcircle parameters
            inRange = [20 50]; resizeFactor = [1 2 5 10];
            
            % noarmalize and rescale images
            objImage = uint8(double(obj.Image)*255/max(double(obj.Image(:))));
            images = arrayfun(@(n)imresize(objImage,n),resizeFactor,'Uni',0);
            
            %---------------------------------------------------------------
            %         Step 1 Find circles using imfindcircles
            %---------------------------------------------------------------
            C(1:numel(resizeFactor)+1,1) = {zeros(0,2)}; R = cell(numel(resizeFactor),1); M = R;
            
            n = 1;
            for m = resizeFactor
                [C{n}, R{n}, M{n}] = imfindcircles(images{n},inRange,'Sensitivity',0.85);
                
                % rescale center and radius
                C{n} = C{n}/m; R{n} = R{n}/m;
                
                n = n+1;
            end
            
            % clear circles with low metric
            centers = cell2mat(C);
            radii = cell2mat(R);
            metric = cell2mat(M);
            
            centers(metric < 0.95,:) = [];
            radii(metric < 0.95) = [];
            metric(metric < 0.95) = [];
            
            % clear overlapping circles
            [m, n] = meshgrid(1:numel(radii), 1:numel(radii));
            X = centers(:,1); Y = centers(:,2);
            overlap = any( ((X(m) - X(n)).^2 + (Y(m) - Y(n)).^2 <= radii(m).^2) & metric(m) > metric(n) ,2);
            
            centers(overlap,:) = [];
            radii(overlap) = [];
            % metric(overlap) = [];
            
            %---------------------------------------------------------------
            %         Step 2 Remove found circles and try again
            %---------------------------------------------------------------
            
            % remove cirlces already found
            Im = images{1};
            [colums, rows] = meshgrid(1:size(Im,2),1:size(Im,1));
            inCirc = @(c,r)(colums-c(1)).^2 + (rows-c(2)).^2 <= (r+2.5).^2;
            for i = 1:length(radii)
                Im(inCirc(centers(i,:),radii(i))) = 0;
            end
            
            A = double(Im); A = A/median(A(A>0));
            B = A; B(Im<25) = 0; % B(I>200) = 255; % potentially do without
            D = imopen(logical(B),strel('disk',5)); % imextendedmax(B,5);

            % find coarse estimate of droplets and filter bad shapes
            % (conjoined droplets)
            imageStats = regionprops(D,A,'BoundingBox','Extent','Image');
            badEstimates = [imageStats.Extent] < 0.75; blank = zeros(size(A)); % start blank guess for droplet locations
            % Note: perfect circles should have extend around 0.78, larger
            % values usually occur for very small droplets

            % filter out conjoined droplets and find a better estimate
            % based on distance maxima
            localRadii = bwdist(~D);
            localMaxima = logical(imregionalmax(localRadii).*bwmorph(imregionalmax(medfilt2(localRadii)),'thicken',2)); % finds maxima and roughly filters local maxima between two droplets
            for i = find(~badEstimates) % remove regions of good estimates from image
                n = imageStats(i);
                Ys = ceil(n.BoundingBox(1)):floor(n.BoundingBox(1)+n.BoundingBox(3));
                Xs = ceil(n.BoundingBox(2)):floor(n.BoundingBox(2)+n.BoundingBox(4));
                localMaxima(Xs,Ys) = zeros(size(n.Image));
                blank(Xs,Ys) = n.Image; % fill up blank guess with good estimates for droplets
            end

            % find located maxima and draw into blank 
            localStats = regionprops(localMaxima,'Centroid','PixelIdxList');
            centroid = vertcat(localStats.Centroid);
            radius = arrayfun(@(x)mean(localRadii(x.PixelIdxList))-3,localStats);

            [colums, rows] = meshgrid(1:size(A,2),1:size(A,1));
            inCirc = @(c,r)(colums-c(1)).^2 + (rows-c(2)).^2 < r.^2;

            for n = 1:length(radius)
                blank(inCirc(centroid(n,:),radius(n))) = 1;
            end
            
            % get sizes of shapes not connected to inlet/outlet
            finalGuess = bwmorph(blank,'thicken',2);
            F = double(A).*finalGuess; F = F/max(F(:)) > exp(-1.5);
            finalStats = regionprops(logical(imclearborder(F)), obj.Image,'WeightedCentroid','Eccentricity','EquivDiameter','EulerNumber','Image');
            finalStats([finalStats.Eccentricity] >= 0.75) = [];
            finalStats([finalStats.EulerNumber] < 1) = [];
            finalStats([finalStats.EquivDiameter]/2 < 2) = [];
            
            %---------------------------------------------------------------
            %         Step 3 Consolidate and report results
            %---------------------------------------------------------------
            warning('on')
            
            centers = [centers; vertcat(finalStats.WeightedCentroid)];
            radii = [radii; vertcat(finalStats.EquivDiameter)/2];

            analysis = struct('for','circles',...
                                    'centers',centers,...
                                    'radii',radii);
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