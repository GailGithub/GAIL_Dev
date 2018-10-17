function x = digitalseq_b2g(s, n)
% function x = digitalseq_b2g(s, n)
%
% Generate the points from a digital sequence in base 2 in gray coded radical
% inverse ordering.
%
% Inputs:
%   s       the number of dimensions per vector
%   n       the number of samples you want
% Outputs:
%   x       an array of sample vectors, size s-by-n
%
% Usage:
%
%   1. Initialize the generator for a point set with generating matrices in the
%   array Cs of dimension s-by-m (see below for the format of the generating
%   matrices array)
%
%     digitalseq_b2g('init0', Cs)
%
%   Valid options for intialization are:
%     init0       the digital sequence as is, the first point is the zero vector
%     init1       the first point is changed into an all ones vector
%     initskip    skip the first point
%
%   2. Generate the next n s-vectors of the sequence, returning an array of
%   dimensions s-by-n:
%
%     P = digitalseq_b2g(s, n)
%
%   3. Ask the state of the current generator with
%
%     digitalseq_b2g('state')
%
%
% Format of generating matrices array:
%
% For a generating matrix in dimension j we define the columns c_1, ..., c_m as
%         [ c_{1,1}  c_{1,2}   ...  c_{1,m}  ]
%         [ c_{2,1}  c_{2,2}   ...  c_{2,m}  ]
%   C_j = [   ...      ...     ...    ...    ] = [ c_1 c_2 ... c_m ]         (*)
%         [ c_{64,1} c_{64,2}  ...  c_{64,m} ]
% This is a binary matrix of dimensions 64-by-m (for standard digital nets these
% matrices are m-by-m, for higher order nets they typicall have the form
% alpha*m-by-m for some integer (interlacing) factor alpha).
% Now represent each column c_k as an integer:
%   c_k = sum_{i=1}^{64} c_{i,k} 2^{i-1},     for k = 1, ..., m.
% I.e., least significant bits are in the top rows.
% These integers are then placed in a row vector to represent C_j like in (*).
% The array of generating matrices is then just an s-by-m array with these
% integer representations as the rows.
% After initialization the radical inverse of these integer column
% representations are taken such that the coordinate x_j can be build by xoring
% these columns and then scaling.
%
%
% Example:
%   load nxmats/nx_b2_m30_s12_Cs.col
%   digitalseq_b2g('init0', nx_b2_m30_s12_Cs)
%   P = digitalseq_b2g(12, 1024);
%   plot(P(5,:), P(12,:), 'r+');
%   axis square;
%   hold on
%   P = digitalseq_b2g(12, 1024);
%   plot(P(5,:), P(12,:), 'b+');
%
% Other nice (actually bad) projections are: 2 versus 5, 4 versus 12, 7
% versus 8, 7 versus 9, 7 versus 10, 8 versus 10, 9 versus 10
%
% (w) 2010, Dirk Nuyens, Department of Computer Science, KU Leuven, Belgium
%     2015 use 64 bit integers and adjusted documentation

persistent k s_max initmode n_max cur recipd Csr maxbit

if ischar(s) && strncmp(s, 'init', 4)
    Cs = n; % when intializing we expect the generating matrices as argument 2 (which is n)
    m = size(Cs, 2);
    s_max = size(Cs, 1);
    n_max = bitshift(uint64(1), m);
    Csr = bitreverse64(uint64(Cs));
    initmode = 0;
    maxbit = 64;
    recipd = pow2(-maxbit);
    k = 0;
    cur = zeros(s_max, 1, 'uint64');
    if strcmp(s, 'init0')
        k = 0;
        initmode = 0;
    elseif strcmp(s, 'init1')
        k = 0;
        initmode = 1;
    elseif strcmp(s, 'initskip')
        initmode = -1;
        k = 1;
    else
        error('I only know about ''init0'', ''init1'' and ''initskip'', what are you talking about?');
    end;
    return;
elseif ischar(s) && strcmp(s, 'state')
    x.index_of_next_point = k;
    x.previous_point_as_integer = cur';
    x.s_max = s_max;
    x.maxbit = maxbit;
    x.n_max = n_max;
    x.Csr = Csr;
    return;
end;

if ((k + n) > n_max) || (s > s_max)
    error('Can only generate %d points in %d dimensions', n_max, s_max);
end;

x = zeros(s, n);

if (k == 0) && (initmode == 0)
    x(:, 1) = 0; si = 2; k = k + 1;
elseif (k == 0) && (initmode == 1)
    x(:, 1) = 1; si = 2; k = k + 1;
else
    si = 1;
end;

for i=si:n
    c = 1;
    while bitget(k, c) == 0
        c = c + 1;
    end;
    cur = bitxor(cur, Csr(1:s_max, c));
    x(:, i) = double(cur(1:s)) * recipd;
    k = k + 1;
end;
