function [freq,varargout] = SAMM_methods(Y,X,nfft,fs,varargin)
%% === example ===
% [freq,eigenvalues,modes] = SAMM_methods(Y,X,nfft,fs) 
% returns the POD modes for all frequencies and all snapshots using SAMM-RR method.
% Y: velocity matrix, which dimension is grids x snapshots;
% X: pressure structure, X(i).pp is the ith sensor matrix, which dimension is snapshots x time delays
% nfft: number of samples for fft calculation
% fs: sampling rate

% Examples:
% [freq,eigenvalues,modes] = SAMM_methods(Y,X,nfft,fs,f_index) 
% f_index: specified frequnecy indices (calculate all frequencies if not specified)

% [freq,eigenvalues,modes] = SAMM_methods(Y,X,nfft,fs,f_index,'SPOD',blocks) 
% blocks: specified number of blocks (use all blocks for SPOD if not specified)
% specify method 'RR' (default) or 'SPOD'

% [freq,RR_eigenvalues,RR_modes,SPOD_eigenvalues,SPOD_modes] = SAMM_methods(Y,X,nfft,fs,f_index,'Both',blocks) 
% This retunrs POD modes for both SAMM-RR and SAMM-SPOD methods

%% === Parameters ===
nsr         = length(X); % number of sensors
[grids,nsp] = size(Y); % grids = resolution; nsp = number of snapshots
nlags       = size(X(1).pp,2); % number of time step in X data
f_center    = nlags/2+0.5; % time step index corresponding to PIV snapshot
vel_std    = std(Y,0,2);
nz          = single(vel_std~=0).*(1:grids)'; % non-zeros vector index
nz(nz == 0) = []; % only keep non-zeros index, increase speed
freq        = 0:fs/nfft:fs/2; % frequency bins
ang_correct = exp(1i*pi*((1:nfft)-1))'; % angle correction for FFT

%% === check the number of inputs ===
method = 'RR'; % default method using RR
blocks = nsp; % default using all snapshots for SAMM-SPOD
f_index = 1:length(freq); % default using all frequencies

if nargin > 4
    f_index = varargin{1}; % specified frequencies
    if nargin == 6
    method  = varargin{2}; % specified SAMM method
    end
    if nargin == 7
    blocks  = varargin{3}; % specified blocks for SAMM-SPOD
    end   
elseif varagin < 4
    error('Insufficient number of inputs');
end

%% === P(f) === 
disp('Calculate G_xx matrix')
p_tmp = zeros(nfft,4772,nsr);
for i = 1:nsp
    for j = 1:nsr
        p_data = X(j).pp(i,f_center-nfft/2:f_center+nfft/2-1); % extract the data
        p_tmp(:,i,j) = sqrt(1/(nfft*fs))*1.59*fft(hamming(nfft).*p_data').*ang_correct; % 1.59 is the correction factor for hamming window
    end
end

pp_fft = zeros(nsr,nsr,nfft);
for ii = 1:nsr
    for jj = 1:ii % only calculate lower triangle
        pp_fft(ii,jj,:) = mean(p_tmp(:,:,ii).*conj(p_tmp(:,:,jj)),2);
        if jj ~= ii 
            pp_fft(jj,ii,:) = conj(pp_fft(ii,jj,:)); % Mapping to upper triangle
        end
    end
end
pp_fft(:,:,2:end-1) = 2*pp_fft(:,:,2:end-1); % two-side psd correction

%% ==== Calculate Rxy matrix ===
disp('Calculate R_xy matrix')
pu = crossR(X,Y,nlags); % cross-correlation between input and output

%% === Calculate Gww and A ===
disp('Calculate G_ww and A')
Gww = zeros(nfft/2+1,nsr);
A = zeros(nfft/2+1,nsr,nsr);
for i = 1:nfft/2+1
    pp_tmp = squeeze(pp_fft(:,:,i));
    [U,S,~] = svd(pp_tmp); %eigendecomposition, SVD speed is comparable to EIG
    Gww(i,:) = diag(S);
    A(i,:,:) = U;
end
if strcmp(method,'RR')
    varargout{1} = Gww;
end
%% === Calculate transfer function Hyw(f) (SAMM-RR) ===
disp('Calculate H_yw')

switch method
    case 'SPOD' 
        H = zeros(grids,nsr,length(f_index));
    case 'RR'                                                                               
        H_yw = zeros(grids,nsr,length(f_index)); 
    otherwise
        H = zeros(grids,nsr,length(f_index));
        H_yw = zeros(grids,nsr,length(f_index));         
end

tmp = zeros(nfft,nsr);
for grid_index = 1:length(nz)
    display(strcat('Processing:',num2str(grid_index/length(nz)*100),'%'))
    for k = 1:nsr % # of sensors    
        tmp(:,k) = fft(hamming(nfft).*reshape(pu(nz(grid_index),k,f_center-nfft/2:f_center+nfft/2-1),[nfft,1]))/fs.*ang_correct; % this is pu_fft, make it temp to save RAM
    end
    tmp(2:end-1,:) = 2*tmp(2:end-1,:); % 2-sided PSD correction 
    for ifreq = 1:length(f_index) % time can be reduced if only calculate specific frequencies frequencies 
        switch method 
            case 'SPOD' 
                H(nz(grid_index),:,ifreq) = (pp_fft(:,:,f_index(ifreq))\tmp(f_index(ifreq),:)')';
            case 'RR'
                H_yw(nz(grid_index),:,ifreq) = (pp_fft(:,:,f_index(ifreq))\tmp(f_index(ifreq),:)')'*squeeze(A(f_index(ifreq),:,:)); % Hyw = HA 
            otherwise
                H(nz(grid_index),:,ifreq) = (pp_fft(:,:,f_index(ifreq))\tmp(f_index(ifreq),:)')';
                H_yw(nz(grid_index),:,ifreq) = (pp_fft(:,:,f_index(ifreq))\tmp(f_index(ifreq),:)')'*squeeze(A(f_index(ifreq),:,:)); % Hyw = HA  
        end
    end
end
if strcmp(method,'RR') || strcmp(method,'Both')
    varargout{2} = H_yw;
end

%% SAMM-SPOD
if strcmp(method,'SPOD') || strcmp(method,'Both')
disp('Calculate SPOD modes')
Phi = zeros(grids,length(f_index));
Gamma = zeros(1,length(f_index));
for ifreq = 1:length(f_index)
    display(strcat('SPOD Processing:',num2str(ifreq/length(f_index)*100),'%'))
    if f_index(ifreq) == 1 || f_index(ifreq) == length(freq)
        Q = squeeze(H(:,1:nsr,ifreq))*squeeze(conj(p_tmp(f_index(ifreq),1:blocks,1:nsr)))';
    else
        Q = sqrt(2)*squeeze(H(:,1:nsr,ifreq))*squeeze(conj(p_tmp(f_index(ifreq),1:blocks,1:nsr)))';
    end
    [Psi,Lambda] = eig(Q'*Q/blocks); % eigen decomposition, method of snapshot
    Lambda       = diag(Lambda);
    [Lambda,idx] = sort(Lambda,'descend');
    Psi          = Psi(:,idx);
    Psi          = Q*Psi*diag(1./sqrt(Lambda)/sqrt(blocks)); % scale to have mode unity length
    Phi(:,ifreq) = Psi(:,1); % Rank-1 mode
    Gamma(ifreq) = Lambda(1); % Rank-1 Eigenvalue
end
end
if strcmp(method,'SPOD')
    varargout{1} = Gamma;
    varargout{2} = Phi;
elseif strcmp(method,'Both')
    varargout{3} = Gamma;
    varargout{4} = Phi;
end    
end

function XY = crossR(X,Y,length_x)
%% calculate cross-covariance of X(input) and Y(output)
% X(sensor).pp(snapshot,time_delay),Y(grid,snapshot)
% length_x = time delay
%%
LX = length(X);
[MY,NY] = size(Y);  % velocity MY = grid,NY = snapshot
XY = single(zeros(MY,LX,length_x));
for i = 1:LX
    temp = single(full(Y)*fliplr(X(i).pp)/NY); 
    XY(:,i,:) = reshape(temp,[MY,1,length_x]); % (grids,sensor,R)
end
end
