% bayesTutorial1
%
% 6/20/98  dhb  Started writing it.
% 7/2/98   dhb  Got through simple linear problem.
% 7/3/98   dhb  Cleaned up, added references.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction
%
% In this tutorial, we will introduce some basic concepts
% related to Bayesian Estimation.  This tutorial focusses
% on a simple example, estimation in the case of a
% well-determined linear problem.  After you study this
% tutorial, you should go onto bayesTutorial2, which
% considers the same ideas in the context of underdeterined
% estimation problems, both linear and non-linear.
% Underdetermined problems are more typical of perceptual
% problems, but it is helpful to introduce the concepts
% with a simple example first.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Problem specification.
%
% Let's suppose there are two pixels in the world,
% so that the scene parameters are described by
% a 2-dimensional column vector.  For this example,
% we'll consider the case where the image is also
% two dimensional thus also described by a 
% 2-dimensional column vector.
%
% In addition, we'll assume that the image is linearly
% related to the scene parameters, but that there is
% some noise.  Thus the rendering equation for this
% problem is:
%
%   imageData = renderingMatrix*sceneParameters + noiseVector
%
% where imageData, sceneParameters, and noiseVector are all
% 2-dimensional column vectors and renderingMatrix is a
% 2 by 2 matrix.

% Let's set a renderingMatrix.  This particular matrix mixes
% the signals from the two pixels, you can think of it as
% a toy version of optical blurring.  Each row of the matrix
% specifies the relation between the sceneParameters and the
% corresponding entry of the imagData.
renderingMatrix = [0.8 0.2 ; 0.2 0.8];

% Let's also set the noise level.  We'll assume independent
% Gaussian noise with zero mean.
noiseMean = [0 0]';
noiseVar = 0.001;
noiseK = [noiseVar 0 ; 0  noiseVar];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LIkelihood.
%
% In the Bayesian formulation, we express the rendering
% equation as a probability distribution, called the
% likelihood.  For any scene parameters, the likelihood is
% the probability, p(imageData | sceneParameters), that
% that we observed imageData.
%
% Let's calculate the likelihood for a particular value
% of imageData for samples of sceneParameters over
% some interesting range.

% Define image data and range of values for scene parameters 
imageData = [0.5 1]';
sceneParameters1Values = 0:0.04:2;
sceneParameters2Values = 0:0.04:2;
nParameters1 = length(sceneParameters1Values);
nParameters2 = length(sceneParameters2Values);
nParameterPairs = nParameters1*nParameters2;

% The meshgrid command creates two matrices
% (here called sceneParameters1 and sceneParameters2)
% whose corresponding entries give all possible pairs
% of the scene parameter values set on each dimension.
% This is a convenient form for plotting but not for calculation.
% So we also create a 2 by nParameterPairs matrix (here
% called sceneParameters) that contains all of pairs strung
% out in its columns.
[sceneParameters1,sceneParameters2] = meshgrid(sceneParameters1Values,sceneParameters2Values);
sceneParameters = [reshape(sceneParameters1,1,nParameterPairs) ; ...
	reshape(sceneParameters2,1,nParameterPairs)];
	
% The likelihood is easily computed given a function
% to compute the probability density function for
% a bivariate normal.  The likelihood is simply
% the probability of a bivariate normal draw 
% being equal to the difference between the
% imageData and the rendered scene parameters,
% where the mean and covariance of the normal
% are given by the mean and covariance of the
% observation noise.  You should be sure you
% understand this, although it may help to
% look at the plot and read the discussion below.
likelihood = ...
	BivariateNormalPdf(imageData(:,ones(1,nParameterPairs))-renderingMatrix*sceneParameters,...
	noiseMean,noiseK);
	
% We'll plot the likelihood as a function of the sceneParameters, given our
% observed imageData.  More generally,the likelihood is a function of both the
% sceneParameters and the imageData, so that for our example, where each of these
% is 2-dimensional, we'd need to plot a function of four variables to see the
% full function.  Note that we use the direct output of meshgrid for the plot,
% and reshape the likelihood into a corresponding form.
figure(1); clf;
mesh(sceneParameters1,sceneParameters2,reshape(likelihood,nParameters2,nParameters1));
xlabel('Scene Parameter 1'); ylabel('Scene Parameter 2'); zlabel('Likelihood');
drawnow;

% The plot is easy to interpret.  When the scene
% parameters are near [0.5 1]', the image data are
% likely to be near [0.5 1]' as well, since given
% our particular renderingMatrix we have
%
%   [0.5 1]' = renderingMatrix*[0.5 1]'
%
% For scene parameters far away from [0.5 1]', observed data of [0.5 1]'
% are not likely. This is because for such case, renderingMatrix*sceneParameters
% is far from the observed [0.5 1]', so you'd need a unusually large noise draw
% to produce an observation of [0.5 1]'.
%
% The spread of the plotted likelihood is due to the noise.
%
% You may want to re-run the code above to see how the likelihood changes if
% you vary imageData or noiseVar.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prior.
%
% In Bayesian estimateion, we express prior constraints about the world
% as a prior distribution.  This specifies how likely we think any value
% of sceneParameters is, before we observe any data.  To get started,
% let's use a Gaussian prior, specified by a mean and covariance matrix.
% We'll compute the prior over the same grid of scene parameters as we
% used for the likelihood.
%
% In this case we choose a rather broad prior.
priorMean = [1 1.5]';
priorK = [10 0 ; 0 10];
prior = BivariateNormalPdf(sceneParameters,priorMean,priorK);
figure(2); clf;
mesh(sceneParameters1,sceneParameters2,reshape(prior,nParameters2,nParameters1));
xlabel('Scene Parameter 1'); ylabel('Scene Parameter 2'); zlabel('Prior');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Posterior.
%
% The posterior is proportional to the product of the likelihood and the prior.
% Generally, we don't need to worry about the normalization constant.
posterior = likelihood .* prior;
figure(3); clf;
mesh(sceneParameters1,sceneParameters2,reshape(posterior,nParameters2,nParameters1));
xlabel('Scene Parameter 1'); ylabel('Scene Parameter 2'); zlabel('Prior');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAP Estimate.
%
% For this problem, it turns out (see comments in function
% NormalPosterior for references) that the posterior is normally
% distributed.  This normality is a consequence of three things.
% First, the problem was linear.  Second, we had normally distributed
% noise. And third, our prior was normal.  For a normal posterior,
% it is reasonable to choose as our best estimate the maximum of
% the posterior, which in this special case is the same as the posterior
% mean.  Choosing the maximum of the posterior is often called MAP
% estimation, where MAP stands for Maximum A Posteriori.
[nil,index] = max(posterior);
mapParameters = sceneParameters(:,index);
fprintf('MAP estimate: (%g,%g)\n',mapParameters(1),mapParameters(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comparison with classical solution.
%
% This problem was not very interesting, as the observation noise
% was low, the rendering matrix was fairly well conditioned,
% and we had fairly weak priors.
%
% Indeed, we could have solved this problem in a non-Bayesian way
% using linear regression.
regParameters = renderingMatrix\imageData;
fprintf('Regression estimate: (%g,%g)\n',regParameters(1),regParameters(2));

% Notice that the regression estimate is quite close to the Bayesian
% estimate.  Indeed, as we'll see when we calculate the Bayesian
% estimate analyitically, most of the difference is due to the 
% quantization of the sceneParameters for the numerical Bayesian
% calculation.  In any case, the close agreement is reassuring.
% For easy problems most methods agree and it probably isn't worth
% fussing much over the differences between them.  We'll come back
% to the difference between the Bayesian approach and simple regression
% when we play with the observation noise below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analytic computation of posterior.
%
% A handy fact to know is that for the case we are considering
% (linear problem, normal noise, normal prior),
% there is a closed form solution for the posterior.  This
% can be derived as part of showing that the posterior is
% in fact normal.  Here we'll just compute the posterior mean.
% The subroutine NormalPosterior encapsulates this calculation
% and also returns the posterior variance.  See the comments
% there to references for the derivations.
W = priorK*renderingMatrix'*inv(renderingMatrix*priorK*renderingMatrix'+noiseK);
b =	priorMean - W*renderingMatrix*priorMean - W*noiseMean;
analyticMapParameters = W*imageData + b;
fprintf('Analytic MAP estimate: (%g,%g)\n',analyticMapParameters(1),analyticMapParameters(2));

% The difference between the explicitly computed and analytic
% estimates is due to the quantization of the underlying scene
% parameter in the case of explicit computation.  You can verify
% this by reducing the sampled range of the scene parameters around
% the maximum and decreasing the spacing betwen samples.  The graphs
% become much less informative in this case, but the numerical estimate
% approaches the analytic one.
%
% The analytic MAP estimate for the linear/normal/normal case is often
% called the discrete Wiener estimate and is often derived from a
% non-Bayesian perspective.  Note that for a normal posterior the MAP
% estimates and the posterior mean are the same and that the mean of
% the posterior (in general) minimizes the expected sum of squares
% estimation error.  It is also worth noting that this estimator is
% often applied to problems that do not satisfy the linear/normal/normal
% conditions for which it is optimal, and that is is not necessarily
% optimal in this case.  For an example, see:
%
%  E. P. Simoncelli and E. H. Adelson, "Noise removal via bayesian
%  wavelet coring," 3rd IEEE Int'l Conf on Image Processing  1-4 (1996)
%
% They explicitly compare an approximately  optimal Bayesian estimator
% to a Wiener estimator for a problem with a non-normal prior.
%
% For general problems, there are not analytic solutions for the
% posterior distribution or even for its maximum/mean.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Effect of prior and noise variance.
%
% With our analytic formula, we can gain a little more insight
% into Bayesian estimate.  Let's reset some of our parameters.
% In particular, we'll make the observations much noisier and
% the examine what happens to the estimate as we vary the 
% prior variance.  We'll use the function NormalPosterior to
% get the posterior max/mean as well as the posterior covariance.
noiseMean = [0 0]';
noiseVar = 0.1;
noiseK = [noiseVar 0 ; 0  noiseVar];
priorMean = [1 1.5]';
priorVars = logspace(2,-4,100);
theEstimates = zeros(2,length(priorVars));
theCovs = cell(length(priorVars));
for i = 1:length(priorVars)
	priorK = [priorVars(i) 0 ; 0 priorVars(i)];
	W = priorK*renderingMatrix'*inv(renderingMatrix*priorK*renderingMatrix'+noiseK);
	b =	priorMean - W*renderingMatrix*priorMean - W*noiseMean;
	[theEstimate,theCov] = NormalPosterior(imageData,renderingMatrix,priorMean,priorK,noiseMean,noiseK);
	theEstimates(:,i) = theEstimate;
	theCovs{i} = theCov;
end
close all
figure(1); clf
subplot(1,2,1);
plot(theEstimates(1,:)',theEstimates(2,:)','g+');
hold on
plot(regParameters(1),regParameters(2),'b*');
plot(priorMean(1),priorMean(2),'r*');
axis([ sceneParameters1Values(1) sceneParameters1Values(length(sceneParameters1Values)) ...
		   sceneParameters2Values(1) sceneParameters2Values(length(sceneParameters2Values)) ]);
xlabel('Scene Parameter 1'); ylabel('Scene Parameter 2');
hold off
subplot(1,2,2);
plot(-log10(priorVars),theEstimates(2,:)','g+');
hold on
plot(-log10(priorVars),theEstimates(2,:)','g');
plot([-log10(noiseVar) -log10(noiseVar)]', ...
	[sceneParameters2Values(1) sceneParameters2Values(length(sceneParameters2Values))]','r');
xlabel('Minus log prior variance'); ylabel('Estimate of parameter 2');
axis([ min(-log10(priorVars)) max(-log10(priorVars)) ...
       sceneParameters2Values(1) sceneParameters2Values(length(sceneParameters2Values)) ]);
drawnow;

% The plots show how the Bayesian estimate changes as
% a function of the strength of the prior.
%
% In the left plot the blue * shows the estimate obtained
% from linear regression.  This estimate ignores both
% prior information and the strength of the observation noise.
% For a broad prior (high variance), the Bayesian 
% estimate agrees with this.  The red * shows the prior
% mean.  As the strength of the prior grows (i.e. as the
% prior variance shrinks), the Bayesian estimate is
% drawn closer and close to the prior mean. 
%
% In the right plot, you can see more clearly how the
% estimate varies as a function of the prior noise.
% Here only the estimate for the second scene parameter
% is plotted as a function of the minus log of the prior 
% variance.  Large prior variances are thus on the left
% of the plot, small prior variances on the right.  The
% red vertical line indicates the variance of the
% observations.
%
% You should think about whether this behavior makes sense
% and whether it is a desirable property for estimation.
% Try increasing the strength of the data by descreasing
% the noise variance.  How does this affect the 
% behavior of the Bayesian estimator?  In the right hand
% plot, can you make intuitive sense of how the location
% of the steep shift in estimate depends on the noise
% variance?

