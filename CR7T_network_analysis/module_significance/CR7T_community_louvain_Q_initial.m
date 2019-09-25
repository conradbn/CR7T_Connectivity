function [M,Q,Qc] = CR7T_community_louvain_Q_initial(W,gamma,M0,Mc_ind,B)

W=double(W);                                % convert to double format
n=length(W);                                % get number of nodes
s=sum(sum(W));                              % get sum of edges

if ~exist('B','var') || isempty(B)
    type_B = 'modularity';
elseif ischar(B)
    type_B = B;
else
    type_B = 0;
    if exist('gamma','var') && ~isempty(gamma)
        warning('Value of gamma is ignored in generalized mode.')
    end
end
if ~exist('gamma','var') || isempty(gamma)
    gamma = 1;
end

if strcmp(type_B,'negative_sym') || strcmp(type_B,'negative_asym')
    W0 = W.*(W>0);                          %positive weights matrix
    s0 = sum(sum(W0));                      %weight of positive links
%     if s0                                   %positive modularity *** ADDED THIS TO AVOID NANs WHEN ONLY NEGATIVE WEIGHTS ***
        B0 = W0-gamma*(sum(W0,2)*sum(W0,1))/s0; 
%     else
%         B0 = 0;
%     end
    W1 =-W.*(W<0);                          %negative weights matrix
    s1 = sum(sum(W1));                      %weight of negative links
    if s1                                   %negative modularity
        B1 = W1-gamma*(sum(W1,2)*sum(W1,1))/s1;
    else
        B1 = 0;
    end
elseif min(min(W))<-1e-10
    err_string = [
        'The input connection matrix contains negative weights.\nSpecify ' ...
        '''negative_sym'' or ''negative_asym'' objective-function types.'];
    error(sprintf(err_string))              %#ok<SPERR>
end
if strcmp(type_B,'potts') && any(any(W ~= logical(W)))
    error('Potts-model Hamiltonian requires a binary W.')
end

if type_B
    switch type_B
        case 'modularity';      B = (W-gamma*(sum(W,2)*sum(W,1))/s)/s;
        case 'potts';           B =  W-gamma*(~W);
        case 'negative_sym';    B = B0/(s0+s1) - B1/(s0+s1);
        case 'negative_asym';   B = B0/s0      - B1/(s0+s1);
        otherwise;              error('Unknown objective function.');
    end
else                            % custom objective function matrix as input
    B = double(B);
    if ~isequal(size(W),size(B))
        error('W and B must have the same size.')
    end
end
if ~exist('M0','var') || isempty(M0)
    M0=1:n;
elseif numel(M0)~=n
    error('M0 must contain n elements.')
end

[~,~,Mb] = unique(M0);
M = Mb;

B = (B+B.')/2;                                          % symmetrize modularity matrix
Hnm=zeros(n,n);                                         % node-to-module degree
for m=1:max(Mb)                                         % loop over modules
    Hnm(:,m)=sum(B(:,Mb==m),2);
end

Q0 = -inf;
Q = sum(B(bsxfun(@eq,M0,M0.')));                        % compute modularity

%% Compute modularity contribution for each module
for ii = 1:numel(Mc_ind)
    nodes = M0 == Mc_ind(ii);
    Bc = B(nodes,nodes);
    Qc(ii,1) = sum(Bc(:));
end





