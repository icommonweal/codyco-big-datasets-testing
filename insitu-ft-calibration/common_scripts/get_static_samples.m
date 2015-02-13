function [ static_samples_len, varargout ] = get_static_samples( varargin )
% Get static samples from a bunch of synchronized samples
% 
% Usage: [static_samples_signal1] = get_static_samples(signal1,derivative1,tol1,..)
if( mod(size(varargin,2),3) ~= 0 )
    fprintf('get_static_samples input error');
    return
end

nr_of_datasets = size(varargin,2)/3;

signal = {};
derivative = {};
tol = [];

nr_of_samples = size(varargin{1},1);

for i = 1:nr_of_datasets
    signal{i} = varargin{3*(i-1)+1};
    derivative{i} = varargin{3*(i-1)+2};
    tol(i) = varargin{3*(i-1)+3};
    assert(size(signal{i},1) == nr_of_samples);
    assert(size(derivative{i},1) == nr_of_samples);
end

is_sample_static = true(nr_of_samples,1);

% Compute samples that violate limits
for i=1:nr_of_datasets
    is_sample_static = and(is_sample_static,all(abs(derivative{i})<tol(i),2));
end

% Compute static samples
static_indices = find(is_sample_static == true);

% Find consecutive sequences of at least N static samples
N = 100;

% Find group of consecutive indices
static_indices_diff=diff(static_indices);
size(static_indices_diff)

b = find([static_indices_diff' inf]>1);
sequence_lenghts = diff([0 b]);

sequence_endpoints_indices = cumsum(sequence_lenghts);

sequence_endpoints = static_indices(sequence_endpoints_indices);

d = sequence_endpoints_indices+1;

sequence_startpoints_indices = [1 d(1:(end-1))];

sequence_startpoints = static_indices(sequence_startpoints_indices);

static_samples_len = [];
varargout = cell(1,nr_of_datasets);
published_static_samples = 1;
for i = 1:size(sequence_lenghts,2)
    if( sequence_lenghts(i) > N )
        static_samples_len(published_static_samples) = sequence_lenghts(i);
        sequence_range = sequence_startpoints(i):sequence_endpoints(i);
        for j = 1:nr_of_datasets
            varargout{j}(published_static_samples,:) = mean(signal{j}(sequence_range,:));
        end
        published_static_samples = published_static_samples+1;
    end
end


end

