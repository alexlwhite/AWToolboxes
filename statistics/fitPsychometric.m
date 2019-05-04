%% function r = fitPsychometric(intensity, responseCorrect, opt)
% Fits a psychometric function of either Gumbel or Weibull form, using the
% Palamedes toolbox (http://www.palamedestoolbox.org/). 
% 
% Inputs: 
% - intensity: a 1xT vector of stimulus intensity values on each of T trials 
%
% - responseCorrect: a 1xT vector of Boolean values (0 or 1) indicating
% whether the subject's response was correct on each of T trials 
%
% - opt: a structure of options for fitting. 
%   Fields of opt can be: 
%
%   * logIntensity: if true, this function will log10 transform the input intensitivty
%   values, and then use the Gumbel psychometric function. otherwise,
%   intensity values are not transformed and the Weibull function is used.
%   Default = true;
%     see http://www.palamedestoolbox.org/weibullandfriends.html
%
%   * fixedSlope: this could be a single number at which to fix the
%   function slope. If NaN, the slope is a free parameter. 
%     Default = NaN;
%  
%   * guessRate: the lower asymptote on the y-axis of p(correct). 
%      Default  = 0.5. 
%   
%   * maxLambda: the maximum value of the lambda parameter, which is 1 minus
%   the upper asymptote. The lapse rate (probability of a total guess not based on the stimulus) is
%   2*lambda. 
%      Default = 0.125;
%
%   * nFitSim: how many repetitions to run in PAL_PFML_GoodnessOfFit. If 0, that function is not run. 
%       Default = 0;
% 
%   * alphaSearchLimits: limits on the initial search grid for the alpha
%      parameter (threshold). Default is set by the range of input intensity values. 
% 
%   * alphaSearchLimits2: limits for the 2nd try if the 1st function fit
%      fails
% 
%   * betaSearchLimits: limits on initial search grid for the beta
%   parameter (slope). Defualt is [0.2 5];
% 
%   * Other parameters describing the number of points in each dimension of
%   the search grids: searchGridNPts, searchGrid2NPts, lambdaSearchNPts,  fitQuailtiyGridNPts
% 
%
% Outputs: 
% r: a structure with fields: 
%   * threshold: 75% correct threshold, in the same scale as the input
%   intensity
% 
%   * fitParams: best fitting parameters [alpha beta gamma lambda] 
% 
%   * lambda: the lambda parameter on its own
% 
%   * opt: the whole opt structure copied in, so you can tell what fitting
%   options were used 
% 
%   * rSqr: simple r^2 for the fit quality 
% 
%   * dev: Deviance (transformed likelihood ratio comparing fit of
%   psychometric function to fit of saturated model), as computed by
%   PAL_PFML_GoodnessOfFit. Is NaN if nFitSim=0;
%   
%   * pdev: proportion of the B Deviance values from simulations that were
%        greater than Deviance value of data. The greater the value of pDev,
%        the better the fit. Is NaN if nFitSim=0;
%        "By somewhat arbitrary convention, researchers agree that the
%        fit is unacceptably poor if pDev is less than 0.05" (Kingdom &
%        Prins, page 73)
%
%
% By Alex L. White, University of Washington 2019 
% 
%  
% 
function r = fitPsychometric(intensity, responseCorrect, opt)

%% Set default options
if nargin<3 && ~exist(opt,'var')
    opt = struct;
end

%whether to log10 transform the intensity values; in which case the
%psychometric function is Gumbel
if ~isfield(opt,'logIntensity')
    opt.logIntensity = true;
end

%Value to fix the slope to. If NaN, slope is free
if ~isfield(opt,'fixedSlope')
    opt.fixedSlope = NaN;
end

%guess rate (lower asymptote)
if ~isfield(opt,'guessRate')
    opt.guessRate = 0.5;
end

%maximum lambda parameter (1-upper asymptote)
if ~isfield(opt,'maxLambda')
    opt.maxLambda = 0.125;
end

%number of simulations for goodness of fit testing
if ~isfield(opt,'nFitSim')
    opt.nFitSim = 0;
end

%other paramters related to inintial search grids
if ~isfield(opt,'alphaSearchLimits')
    rng = [min(intensity) max(intensity)];
    opt.alphaSearchLimits = rng + [-1 1]*0.05*diff(rng);
    opt.alphaSearchLimits(opt.alphaSearchLimits<=0) = 0.001*max(intensity);
end
if ~isfield(opt,'alphaSearchLimits2')
    opt.alphaSearchLimits2 = opt.alphaSearchLimits+[-1 1]*0.05*diff(opt.alphaSearchLimits);
    opt.alphaSearchLimits2(opt.alphaSearchLimits2<=0) = 0.001*max(intensity);
end
if ~isfield(opt,'betaSearchLimits')
    opt.betaSearchLimits = [0.2 5];
end

if ~isfield(opt,'searchGridNPts')
    opt.searchGridNPts = 200;
end
if ~isfield(opt,'searchGrid2NPts')
    opt.searchGrid2NPts = 400;
end
if ~isfield(opt,'lambdaSearchNPts')
    opt.lambdaSearchNPts = 10;
end
if ~isfield(opt,'fitQuailtiyGridNPts')
    opt.fitQuailtiyGridNPts = 100;
end

%% transform data? 
if opt.logIntensity
    intensity = log10(intensity);
    psychometric = @PAL_Gumbel;
    
    opt.alphaSearchLimits = log10(opt.alphaSearchLimits);
    opt.alphaSearchLimits2 = log10(opt.alphaSearchLimits2);
else
    psychometric = @PAL_Weibull;
end

%% set up fit options 
lapseLimits = [0 opt.maxLambda];

%structure defining grid to search for initial values:
searchGrid.alpha = linspace(opt.alphaSearchLimits(1), opt.alphaSearchLimits(2),opt.searchGridNPts);
searchGrid.beta  = linspace(opt.betaSearchLimits(1),  opt.betaSearchLimits(2), opt.searchGridNPts);
searchGrid.gamma = opt.guessRate; 
searchGrid.lambda = linspace(0,opt.maxLambda, opt.lambdaSearchNPts); 

%2nd try search grid - different starting levels
searchGrid2.alpha = linspace(opt.alphaSearchLimits2(1), opt.alphaSearchLimits2(2),opt.searchGrid2NPts);
searchGrid2.beta  = linspace(opt.betaSearchLimits(1),   opt.betaSearchLimits(2), opt.searchGrid2NPts);
searchGrid2.gamma = opt.guessRate;  
searchGrid2.lambda = linspace(0,opt.maxLambda, opt.lambdaSearchNPts); 

%search grid for goodness-of-fit: 
searchGridGoF.alpha = linspace(opt.alphaSearchLimits(1), opt.alphaSearchLimits(2),opt.fitQuailtiyGridNPts);
searchGridGoF.beta = linspace(opt.betaSearchLimits(1),   opt.betaSearchLimits(2), opt.fitQuailtiyGridNPts);
searchGridGoF.gamma = opt.guessRate; %guess rate (fixed)
searchGridGoF.lambda = linspace(0,opt.maxLambda, opt.lambdaSearchNPts); 


if ~isnan(opt.fixedSlope)
    searchGrid.beta = opt.fixedSlope;
    searchGrid2.beta = opt.fixedSlope;
    searchGridGoF.beta = opt.fixedSlope;
end

%        [threshold slope guess-rate lapse-rate]
freeParams = [1 isnan(opt.fixedSlope) 0 1]; %2nd parameter is slope, which is free if input fixedSlope is NaN


%% prepare the data: compute number of trials and p(correct) at each intensitiy level 
% don't do any binning except when two trials have exactly the same
% intensity level 

uX = unique(intensity);
nBin = length(uX);
nc = zeros(1,nBin);
nt = zeros(1,nBin);
for i=1:nBin
    binTs = find(intensity==uX(i));
    nc(i) = sum(responseCorrect(binTs));
    nt(i) = length(binTs);
end

pc = nc./nt; %p(correct)

%% Fit 
[fitParams, ~, exitflag] = PAL_PFML_Fit(uX, nc, nt, searchGrid, freeParams, psychometric,'lapseLimits',lapseLimits);

r.firstTryWorked = exitflag;

%if fit didn't converge the first time, try again
if exitflag == 0 
    [fitParams, ~, exitflag] = PAL_PFML_Fit(uX, nc, nt, searchGrid2, freeParams, psychometric,'lapseLimits',lapseLimits);
    r.secondTryWorked = exitflag;
end

%% compute 75% correct threshold
thresh = psychometric(fitParams,0.75,'Inverse');

if opt.logIntensity
    thresh = 10^thresh;
end

%% goodness of fit
%r^2
resid=pc-psychometric(fitParams,uX);
SSTot=sum((pc-mean(pc)).^2);
SSE=sum(resid.^2);
rSqr=1-SSE/SSTot;

%pdev
if opt.nFitSim>0
    [dev, pdev] = PAL_PFML_GoodnessOfFit(uX, nc, nt, fitParams, freeParams, opt.nFitSim, psychometric,'searchGrid', searchGridGoF,'lapseLimits',lapseLimits);
else
    dev = NaN; pdev = NaN;
end
%this function compares our "target" model (from fitting in the previous step) to a "saturated model" in
%which accuracy at each stimulus level is a free parameter. Target model is
%nested inside the saturated model. This function simulates an observer
%using the target model, then fits with the assumptions of the target model
%(certain nPFBins shape with some free parameters), and fits with the less
%restrictive assumptions of the saturated model. For each fit, a
%likelihood, and for each simulation, a likelihood ratio. Likelihood for
%saturated model will always be greater, can always get closer.
%If likelihood ratio from experimental data is often lower than what comes
%from simulations, then some assumptions made by target model are bad.
%"dev": Deviance (transformed likelihood ratio comparing fit of psychometric function to fit of saturated model)
%"pdev": proportion of the B Deviance values from simulations that were
%        greater than Deviance value of data. The greater the value of pDev,
%        the better the fit.
%        "By somewhat arbitrary convention, researchers agree that the
%        fit is unacceptably poor if pDev is less than 0.05" (Kingdom &
%        Prins, page 73)


%% save results 
r.threshold = thresh;
r.fitParams = fitParams;
r.lambda = fitParams(4);
r.opt = opt;
r.rSqr = rSqr;
r.dev = dev;
r.pdev = pdev;

