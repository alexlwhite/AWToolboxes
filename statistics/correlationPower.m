%% Correlation Power Analysis
% Based on an R script by  Guillaume A Rousselet: https://github.com/GRousselet/blog/blob/master/corrpower/docs/corr_power.md#power-as-a-function-of-sample-size-and-rho-homoscedastic-case
% https://garstats.wordpress.com/2018/06/29/corrpower/ 
% 
% Inputs: 
% -  rho: true assumed "rho", Pearson's linear correlation coefficient between two homoscedastic
%    variables.
% - goalPower: power you want. That is, the probability of getting p<0.05
%    in your linear correlation, assuming rho. Default = 0.8;
% - Ns: a vector of sample sizes you want to simulate. Deault: 10:5:200. 
% - nsim: how many simulations to run for each sample size. Default: 10000.
% - doPlot: boolean, whether to make a plot. Default = true;

% Outputs: 
% - minN: the minimum sample size that achieves power>=goalPower. 

function minN = correlationPower(rho, goalPower, Ns, nsim, doPlot)

if nargin<2
    goalPower = 0.8; 
end

if nargin<3
    Ns = 10:5:200;
end
nN = length(Ns);

if nargin<4
    nsim = 10000;
end

if nargin<5
    doPlot = true;
end



corr_sig = zeros(1, nN);

for ni = 1:nN
    n = Ns(ni);

    % Define standard deviations and means
    sd_x = ones(1, n);
    sd_y = ones(1, n);
    mu_x = zeros(1, n);
    mu_y = zeros(1, n);

    for iter = 1:nsim

        % Generate data using Y conditional on X
        % https://www.r-bloggers.com/simulating-from-the-bivariate-normal-distribution-in-r/
        % Proof: http://athenasc.com/Bivariate-Normal.pdf
        x = normrnd(mu_x, sd_x);
        y = normrnd(mu_y + (sd_y ./ sd_x) .* rho .* (x - mu_x), sqrt(1 - rho^2) .* sd_y.^2);
        [~, pval] = corr(x', y');

        % Update power count
        if pval <= 0.05
            corr_sig(ni) = corr_sig(ni) + 1;
        end
    end
end

powers = corr_sig/nsim;

enough = find(powers>=goalPower);
minN = min(Ns((enough))); 

if isempty(minN)
    fprintf(1,'\nYour set of sample sizes didn''t go big enough!!\n');
end

if doPlot

figure; hold on;
plot(Ns, powers,'.-');
if ~isempty(minN)   
    plot([minN minN], [0 goalPower], '-');
    ttl = sprintf('N=%i needed for %i%% power to detect r=%.2f', minN, round(100*goalPower), rho);
else
    ttl = sprintf('Need bigger Ns for %i%% power to detect r=%.2f', round(100*goalPower), rho);
end
xlabel('Sample size'); ylabel('Proportion p<=0.05');

title(ttl);
xlim([0 max(Ns)*1.1]);
ylim([0 1]);

end

