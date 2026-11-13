function calculate_anomalies(inputFile, variableName, coefficientFile, varargin)
%CALCULATE_ANOMALIES Fit grid-cell time series and save compact coefficients.
%   CALCULATE_ANOMALIES(INPUTFILE, VARIABLENAME, COEFFICIENTFILE) fits an
%   intercept, linear trend, and seasonal harmonics independently at each
%   lon/lat/depth grid cell, reading the input NetCDF in spatial tiles.
%   The output NetCDF contains the fitted coefficients in variable "beta";
%   it does not contain the full anomaly time series.
%
%   The input variable must have exactly three spatial dimensions and one
%   time dimension. Periods and the time coordinate must use the same units.
%   The default seasonal periods assume that time is measured in days.
%
%   Name-value options:
%     'TileSize' - tile size along the three spatial dimensions (default
%                  [8 8 2]).
%     'Periods'  - seasonal periods (default [365.2425 365.2425/2]).
%     'MinValid' - minimum finite observations per grid cell (default is
%                  number of fitted coefficients plus one).
%
%   The beta coefficient order is intercept, linear trend, then sine and
%   cosine coefficients for each period in 'Periods'. The trend uses the
%   input time units, with time centered at the output attribute
%   'time_reference'. To evaluate the fitted baseline at time T, use
%   [1, T-time_reference, sin(2*pi*(T-time_reference)./periods), ...
%   cos(2*pi*(T-time_reference)./periods)] * beta.

if nargin < 2
    error('calculate_anomalies:NotEnoughInputs', ...
        'Input file and variable name are required.');
end
if ~isfile(inputFile)
    error('calculate_anomalies:InputNotFound', ...
        'Input NetCDF file does not exist: %s', inputFile);
end
if nargin < 3 || isempty(coefficientFile)
    [folder, name] = fileparts(inputFile);
    coefficientFile = fullfile(folder, [name '_anomaly_coefficients.nc']);
end
if isfile(coefficientFile)
    error('calculate_anomalies:OutputExists', ...
        'Output file already exists; choose another path: %s', coefficientFile);
end

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'TileSize', [8 8 2], ...
    @(x) isnumeric(x) && isreal(x) && numel(x) == 3 && ...
    all(isfinite(x)) && all(x >= 1) && all(x == fix(x)));
addParameter(parser, 'Periods', [365.2425 365.2425/2], ...
    @(x) isnumeric(x) && isreal(x) && isvector(x) && ...
    all(isfinite(x)) && all(x > 0));
addParameter(parser, 'MinValid', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 1 && x == fix(x)));
parse(parser, varargin{:});
tileSize = double(parser.Results.TileSize(:)');
periods = double(parser.Results.Periods(:)');

sourceInfo = ncinfo(inputFile);
variableIndex = find(strcmp({sourceInfo.Variables.Name}, variableName), 1);
if isempty(variableIndex)
    error('calculate_anomalies:VariableNotFound', ...
        'Variable "%s" was not found in %s.', variableName, inputFile);
end
variableInfo = sourceInfo.Variables(variableIndex);
dimensions = variableInfo.Dimensions;
dimensionNames = {dimensions.Name};
timeDimIndex = find(strcmpi(dimensionNames, 'time'), 1);
if isempty(timeDimIndex) || numel(dimensions) ~= 4
    error('calculate_anomalies:UnsupportedDimensions', ...
        'The selected variable must have exactly three spatial dimensions and one dimension named time.');
end
spatialDimIndices = setdiff(1:numel(dimensions), timeDimIndex, 'stable');
spatialDimensions = dimensions(spatialDimIndices);
spatialSizes = double([spatialDimensions.Length]);
numTimes = double(dimensions(timeDimIndex).Length);

if any(~isfinite(spatialSizes)) || any(spatialSizes < 1) || numTimes < 1
    error('calculate_anomalies:InvalidDimensions', ...
        'The selected variable has an empty or unlimited spatial/time dimension.');
end
tileSize = min(tileSize, spatialSizes);
timeVariable = find_coordinate_variable(sourceInfo, dimensions(timeDimIndex).Name);
if isempty(timeVariable)
    error('calculate_anomalies:TimeCoordinateNotFound', ...
        'Could not find a one-dimensional coordinate for the time dimension.');
end
time = double(ncread(inputFile, timeVariable.Name));
time = time(:);
if numel(time) ~= numTimes || any(~isfinite(time)) || any(diff(time) <= 0)
    error('calculate_anomalies:InvalidTime', ...
        'Time coordinates must be finite, strictly increasing, and match the time dimension.');
end

numCoefficients = 2 + 2*numel(periods);
minValid = parser.Results.MinValid;
if isempty(minValid)
    minValid = numCoefficients + 1;
end
if minValid < numCoefficients
    error('calculate_anomalies:MinValidTooSmall', ...
        'MinValid must be at least the number of fitted coefficients (%d).', ...
        numCoefficients);
end
if numTimes < minValid
    error('calculate_anomalies:InsufficientTimeSteps', ...
        'The file has %d time steps, fewer than MinValid (%d).', ...
        numTimes, minValid);
end

timeReference = mean(time);
tau = time - timeReference;
design = [ones(numTimes, 1), tau];
for p = 1:numel(periods)
    phase = 2*pi*tau/periods(p);
    design = [design, sin(phase), cos(phase)]; %#ok<AGROW>
end
if rank(design) < numCoefficients
    error('calculate_anomalies:RankDeficientDesign', ...
        'The selected time values cannot support the requested trend and seasonal fit.');
end

% Create the output variable using the source spatial dimension names/order.
if ~isempty(fileparts(coefficientFile)) && ...
        ~isfolder(fileparts(coefficientFile))
    mkdir(fileparts(coefficientFile));
end
outputDimensions = cell(1, 2*numel(spatialDimensions) + 2);
for k = 1:numel(spatialDimensions)
    outputDimensions{2*k-1} = spatialDimensions(k).Name;
    outputDimensions{2*k} = spatialSizes(k);
end
outputDimensions(end-1:end) = {'coefficient', numCoefficients};
nccreate(coefficientFile, 'beta', 'Dimensions', outputDimensions, ...
    'Datatype', 'single', 'FillValue', single(NaN), 'DeflateLevel', 4);
nccreate(coefficientFile, 'coefficient', ...
    'Dimensions', {'coefficient', numCoefficients}, 'Datatype', 'int32');
ncwrite(coefficientFile, 'coefficient', int32(1:numCoefficients));
ncwriteatt(coefficientFile, 'beta', 'long_name', ...
    'Grid-cell least-squares fit coefficients');
ncwriteatt(coefficientFile, 'beta', 'time_reference', timeReference);
ncwriteatt(coefficientFile, 'beta', 'seasonal_periods', periods);
ncwriteatt(coefficientFile, 'beta', 'minimum_valid_observations', int32(minValid));
ncwriteatt(coefficientFile, 'beta', 'coefficient_order', ...
    'intercept, linear_trend, then sine/cosine pairs for each seasonal_period');
ncwriteatt(coefficientFile, '/', 'source_file', inputFile);
ncwriteatt(coefficientFile, '/', 'source_variable', variableName);
ncwriteatt(coefficientFile, '/', 'processing_status', 'incomplete');

% Copy small coordinate variables so the coefficient grid remains georeferenced.
for k = 1:numel(spatialDimensions)
    coordinate = find_coordinate_variable(sourceInfo, spatialDimensions(k).Name);
    if ~isempty(coordinate)
        copy_coordinate(inputFile, coefficientFile, coordinate, ...
            spatialDimensions(k).Name, spatialSizes(k));
    end
end
copy_coordinate(inputFile, coefficientFile, timeVariable, ...
    dimensions(timeDimIndex).Name, numTimes);

fillValues = get_missing_values(variableInfo);
starts = ones(1, numel(dimensions));
counts = double([dimensions.Length]);

for i1 = 1:tileSize(1):spatialSizes(1)
    count1 = min(tileSize(1), spatialSizes(1) - i1 + 1);
    for i2 = 1:tileSize(2):spatialSizes(2)
        count2 = min(tileSize(2), spatialSizes(2) - i2 + 1);
        for i3 = 1:tileSize(3):spatialSizes(3)
            count3 = min(tileSize(3), spatialSizes(3) - i3 + 1);
            tileStarts = [i1 i2 i3];
            tileCounts = [count1 count2 count3];
            for k = 1:3
                starts(spatialDimIndices(k)) = tileStarts(k);
                counts(spatialDimIndices(k)) = tileCounts(k);
            end
            starts(timeDimIndex) = 1;
            counts(timeDimIndex) = numTimes;

            values = ncread(inputFile, variableName, starts, counts);
            values = reshape(values, counts);
            values = permute(values, [spatialDimIndices timeDimIndex]);
            values = reshape(double(values), [], numTimes);
            for k = 1:numel(fillValues)
                values(values == fillValues(k)) = NaN;
            end
            values(~isfinite(values)) = NaN;

            beta = nan(size(values, 1), numCoefficients);
            complete = all(isfinite(values), 2);
            if any(complete)
                beta(complete, :) = (design \ values(complete, :)')';
            end
            incompleteRows = find(~complete);
            for row = incompleteRows'
                valid = isfinite(values(row, :));
                if nnz(valid) >= minValid && rank(design(valid, :)) == numCoefficients
                    beta(row, :) = (design(valid, :) \ values(row, valid)')';
                end
            end

            beta = reshape(single(beta), [tileCounts numCoefficients]);
            outputStarts = [tileStarts 1];
            ncwrite(coefficientFile, 'beta', beta, outputStarts);
        end
    end
end

ncwriteatt(coefficientFile, '/', 'processing_status', 'complete');
end

function coordinate = find_coordinate_variable(fileInfo, dimensionName)
coordinate = [];
for k = 1:numel(fileInfo.Variables)
    candidate = fileInfo.Variables(k);
    if isscalar(candidate.Dimensions) && ...
            strcmp(candidate.Dimensions(1).Name, dimensionName)
        coordinate = candidate;
        return
    end
end
end

function copy_coordinate(inputFile, outputFile, coordinate, dimensionName, dimensionLength)
if ~isnumeric(ncread(inputFile, coordinate.Name))
    return
end
nccreate(outputFile, coordinate.Name, ...
    'Dimensions', {dimensionName, dimensionLength}, 'Datatype', 'double');
ncwrite(outputFile, coordinate.Name, double(ncread(inputFile, coordinate.Name)));
for k = 1:numel(coordinate.Attributes)
    attribute = coordinate.Attributes(k);
    if any(strcmp(attribute.Name, {'units', 'long_name', 'standard_name', 'axis'}))
        ncwriteatt(outputFile, coordinate.Name, attribute.Name, attribute.Value);
    end
end
end

function fillValues = get_missing_values(variableInfo)
fillValues = [];
for k = 1:numel(variableInfo.Attributes)
    attribute = variableInfo.Attributes(k);
    if any(strcmp(attribute.Name, {'_FillValue', 'missing_value'})) && ...
            isnumeric(attribute.Value)
        fillValues = [fillValues; double(attribute.Value(:))]; %#ok<AGROW>
    end
end
end
