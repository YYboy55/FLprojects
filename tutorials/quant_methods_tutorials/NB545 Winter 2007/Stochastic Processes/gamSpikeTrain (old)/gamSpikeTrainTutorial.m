%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GamSpikeTrainTutorial%% This exercise illustrates a few simple properties of spike intervals.% It introduces a few facts about renewal processes, beginning with the% poisson point process and extending intuitions to gamma processes.% It shows how regular spikes would be if they were to represent% waiting times to the Nth arrival in a Poisson point process.%% This tutorial was created by M. N. Shadlen 4/22/98  and was modified% by Greg Horwitz (for CSH '00). % 1/24/03 MNS updated for Neubeh545%  % Dependencies: % 	plot1ras: written by Mike Shadlen%	play_spikes: written by Geoff Boynton%   plot2ras: written by Greg Horwitz%%%% Notes to me... Things to do:% Add Markov chain?% Preliminary matlab stuff.% Make sure the utilities are in your pathaddpath('/Users/mike/Documents/Matlab_mfiles/quant_methods_tutorials/Stochastic Processes/gamSpikeTrain')% What is a point process?   It is a list of events that are% indistinguishable from one another, like spikes or radioactive decay, or% the onset of an alarm,  or a grain of silver on a piece of film. We can% think of each *event* as having a time or a location.   That means we can% describe the process by a list of times, sometimes called arrival times.% We would characterize a point process by a list of times, or intervals,% or  the number of events in some epoch.  % What is a renewal process? It is a point process whose intervals are% drawn randomly from a common distribution.  Each interval is independent% from the others.  So we say that the intervals are iid (independent and% identically distributed).  They are described by a random variable. %% Example 1. Here is an example of a renewal process.% We're going to draw interevent intervals from a lognormal% distribution.  A lognormal distribution look like this:x = [-5:.1:50];y = lognpdf(x,2.5,1);y = lognpdf(x,1.88,1.5);plot(x,y);title('The Lognormal Distribution');xlabel('Interval duration (ms)');ylabel('Probability');% Notice that the lognormal density has non-zero probability for only% positive interval durations.  This makes sense because the time% between two events can't possibly be negative -- lengths (whether in % time or space) are always non-negative.% (As aside: if X has a lognormal distribution then log(X) has a normal% (Gaussian) distribution.  That's where the name, lognormal, comes from.)% Let's imagine that a neuron fires spikes according to a renewal process. % Further, let's imagine that the interspike intervals from this neuron % have a lognormal distribution.  We can simulate spikes trains by drawing% a bunch of interspike intervals at random from the lognormal distribution, % and calling each event a spike.  We'll repeat this 100 times to simulate% 100 trials.ntrials = 100intervals = lognrnd(1.88,1.5,100,ntrials);% Here's a histogram of the interspike intervals:figure(1)[n,bins] = hist(intervals(:),x);n = n/(sum(n)*.1);bar(x(1:end-1),n(1:end-1));set(gca,'XLim',[min(x) max(x)]);% (A technical point: we cut off the last bin in the histogram % (the line: "bar(x(1:end-1),n(1:end-1)" above) because the way% the 'hist' function works the rightmost bin contains all of % the interspike intervals that fall into *or exceed* that bin.% There are a bunch of interspike intervals off to the right of% the plot which would have been stuck into the rightmost bin% had we not done this.  The relatively high number of interspike% intervals that are very large is due to the fact that the % lognormal distribution has a very long tail -- a characteristic % that we'll come back to soon.)% As expected, the histogram looks a lot like the lognormal % distribution that the interspike intervals were drawn from.  % Let's superimpose the lognormal density to verify the similarity.hold on;plot(x,y,'g-');hold off;% We can use the 'cumsum' (cumulative sum) function to add up the% interspike intervals and thus find the arrival time of each spike.arrivals = cumsum(intervals);figure(1)for i=1:ntrials	plot1ras(arrivals(:,i),i); hold on;endhold off;xlabel('Time (ms)');ylabel('Trial number');% Notice that the arrivals look pretty irregular.% Let's listen to this irregular point process by playing the events % as clicks in the computer speaker.  This is similar to what an% electrophysiologist hears over the audio monitor when recording % from a cell in vivo.  We'll play the spikes from the first five trials% and wait 1 second between trials.for i = 1:5	play_spikes(arrivals(:,i));	pause(1)end% Now that we've generated some simulated spike trains we can do % some statistical analyses of our simulated data.% First we'll calculate the mean and standard deviation of the % interspike intervals.mean_int = mean(intervals(:))var_int = var(intervals(:))% The ratio of standard deviation to mean is called the % 'coefficient of variation'.  CV_int = sqrt(var_int)/mean_int% Now we compute the theoretical value of the CV (with enough simulated% spikes the empirical CV will tend to this value).[m,v] = lognstat(1,2);sqrt(v)/m% (An aside: The CV for the lognormal(2,1) distribution is '1.3108'.% This is larger than '1.0' which is the CV of the exponential distribution.  % "So?", you ask, "why compare to the exponential distribution?"% Because a renewal process, with exponentially  distributed interevent % intervals, is the well-known Poisson process which has all kinds of% interesting properties.  The differences in CV for these two interevent% distributions has consequences for the statistics of spike counts% as we'll see in a moment.)% Neurophysiologists often calculate statistics from spike count data.% In this next bit we're going to calculate spike count statististics% from our simulated spike trains.% We generated 100 spikes for each of the trials in our simulation.  % Because the spike times are random, this means that each trial% ends at a different time.  We will now calculate spike count statistics% during the first 500 ms of each trial.  This is a brief enough epoch% that even the shortest trial exceeds it.epochDur = 500;if (min(max(arrivals)) < epochDur);	error('Not enough spikes!');endcounts = sum(arrivals < epochDur);% Here is the frequency histogram of spike counts.figure(1)hist(counts,[min(counts):2:max(counts)]);xlabel('Number of spikes in epoch')ylabel('Number of trials');% We can now calculate the mean and variance of the distribution % of spike counts (just as we previously did for the distribution% of interspike intervals).mean_counts = mean(counts)var_counts = var(counts)% For any renewal process there's a remarkable asymptotic relationship% between the mean and variance of the interspike interval distribution% and the mean and variance of the distribution of spike counts.% First, the mean of the spike count is (asymptotically) equal to the % duration of the counting window divided by the mean of the interval% distribution.epochDur/mean_intmean_counts% Second, the variance of the spike count is (asymptotically) equal% to the duration of the counting window times the variance of the % interval distribution, divided by the mean of the interval distribution% cubed.  (Note: this converges a lot more slowly than the formula for the% mean, above).epochDur*var_int/mean_int^3var_counts% Another way of saying this is that, for a renewal process, the square% of the CV of the interval distribution is asymptotically equal to % the Fano factor (variance to mean ratio) of the counts.CV_int^2var_counts/mean_counts% Just before we leave renewal processes in general it's worth pointing % that as the length of time over which the process is obsevered becomes% very long (with respect the mean interarrival time) something cool % happens to the distribution of counts within that window: it becomes% more and more Gaussian.  This is a result of the Central Limit theorem.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Example 2.  Generating a Poisson process.% In this exercise we will use the same simulation technique to % generate spike trains from a few different hypothetical neurons.% The first is a Poisson neuron that fires spikes according to % a Poisson process -- a special type of renewal process.  Then,% we'll look at neurons that fire according to a renewal process% with interarrival times that have a gamma distribution.  Such% a neuron can be thought of as a perfect integrator of Poisson% inputs, as will be explained shortly.%% A Poisson process is a renewal process whose interevent times% have an exponential distribution.  To simulate the spiking of% a Poisson neuron, we are going to draw random numbers from % an exponential distribution and place spikes according to these% random numbers.%% (An aside: Note that this is a different approach to simulating% a Poisson process than is used in the "poissonTutorial.m".  This% is actually a better technique for many purposes.  Importantly,% we are not making the simplifying ssumption that time comes in% discrete steps, as was done in "poissonTutorial.m".)ntrials = 100;spike_rate = 50;                 % spikes/secmeanint = 1000 / spike_rate;     % mean interspike interval (ms)intervals = exprnd(meanint,100,ntrials);% Here's the histogram of the interspike intervals that we've % drawn.  To no one's great surprise, it looks pretty exponential.[n,bins] = hist(intervals(:),x);n = n/(sum(n)*.1);bar(x(1:end-1),n(1:end-1));set(gca,'XLim',[min(x) max(x)]);xlabel('Interval duration (msec)')ylabel('Probability')% Now we add up the interspike intervals to find the arrival time% of each spike.  We're also going to play the spikes over the % speaker and draw the raster to the screen.arrivals = cumsum(intervals);figure(1)for i=1:ntrials	plot1ras(arrivals(:,i),i); hold on;endhold off;xlabel('Time (ms)');ylabel('Trial number');for i = 1:5	play_spikes(arrivals(:,i));	pause(1)end% Notice that the spike rasters look qualitatively similar % to the neuron whose interspike intervals had a lognormal% distribution.  The spike trains sound about the same over the % speaker too.  This is because our eyes and ears aren't so good% at discriminating between these two kinds of random processes.% The statistics of these spike trains are rather different from % those we simulated at the begining of this tutorial.% To illustrate this point, we'll now look at the statistics of the % spike counts across trials.epochDur = 500;if (min(max(arrivals)) < epochDur);	error('Not enough spikes!');endpoisscounts = sum(arrivals < epochDur);% Here is the frequency histogram of spike counts.minbin = min([poisscounts,counts]);maxbin = max([poisscounts,counts]);binwidth = 1;figure(1)subplot(2,1,1)[n,x] = hist(poisscounts,[minbin:binwidth:maxbin]);bar(x,n/sum(n));ylabel('Proportion of trials');title('Poisson neuron','FontSize',14);% For comparison, we'll also display the frequency histogram of% spike counts for the neuron we simulated earlier -- the one% whose interspike intervals come from the lognormal distribution.subplot(2,1,2)[n,x] = hist(counts,[minbin:binwidth:maxbin]);bar(x,n/sum(n));xlabel('Number of spikes in epoch')ylabel('Proportion of trials');title('non-Poisson neuron','FontSize',14);% Notice that the number of spikes fired by the Poisson neuron,% although random from trial to trial, is much more consistent% that the number of spikes fired by the neuron from the begining% of the tutorial.  This is a consequence of the fact that the % coefficient of variation of the lognormal distribution (with the % parameters that we used) is bigger than the coefficient of % variation of the exponential distribution.% OK, so what's the big deal about the Poisson process? % What's so great about having exponentially distributed % interevent times?  % The Poisson process has a number of interesting qualities% two of which we'll mention now.  The first is that % the number of events (spikes in our case) arriving during% an epoch of fixed duration has a distribution that can be% derived analytically: the probability that the neuron fires % 'k' spikes during an interval of length 't' is:%%                    (lambda*t)^k e^(-lambda*t)% Prob (X(t) = k) = --------------------------%                                k!%% To illustrate this, we'll superimpose these analytically% derived probabilities over our histogram.subplot(2,1,1);hold on;x = [minbin:binwidth:maxbin];plot(x,binwidth*poisspdf(x,spike_rate*(epochDur/1000)),'m*');% Notice that this theoretical distribution matches the spike counts% from the Poisson neuron pretty well, but it doesn't match the% spikes counts from the other neuron at all.  Sadly, the spike% counts from the non-Poisson neuron don't have a distribution with% a nice analytical solution, at least as far as I know.% You may have noticed that that superimposed Poisson distribution % looks almost like a normal distribution.  The fact is that the % Poisson distribution converges to a normal distribution as the% duration of the counting window (or the rate of events) goes to % infinity.  So does the distribution of spike coutns from the other% neuron.  As was mentioned briefly earlier, the distribution of% counts from any renewal process converges to a normal distribution % as the number of counts becomes large.% The Poisson process has a remarkable property that falls easily% of the exponentially distributed intervals: the probabillity% of firing a spike does not depend on when the last spike% was fired.  Because of this property, the Poisson process% is called "memoryless".% To illustrate this we'll first plot the exponential % density function.figure(1); clf;t = [0:.1:200];y = exppdf(t,meanint);plot(t,y);xlabel('Interval duration');ylabel('Probability');% Looking at the probability distribution you might think that  if a spike% just occurred another one is probably going to happen pretty soon.  After% all, the greatest probability density is for very short intervals (at the% left of the graph).  The fact is that the probability of a spike% happening doesn't  depend at all on when the last spike happened.  They% are  completely independent, for the Poisson process.%% To see this let's consider a conditional probability: the  probability of% firing a spike at time 't' given that we haven't fired a spike from time% 0 to 't'.  Assume that a spike was fired at time 0 exactly.  This first% probability is just the  exponential probability density (it is the% probability  of observing an interspike interval 't' ms long).  The% second probability is one minus the cumulative exponential density from 0% to 't' (the probability that a spike hasn't occurred between time 0 and% 't').   To find the conditional probability, we divide one by the other.hazard = exppdf(t,meanint)./(1-expcdf(t,meanint));clf;plot(t, hazard)xlabel('Time (ms)')ylabel('Prob spike | Prev spike at t=0')set(gca,'Ylim',[0 .1]);% The so-called "hazard function" is totally flat.  This means that % the probability of firing a spike in the next small time increment% is the same irrespective of how much time has gone by since the % previous spike.  The product of the hazard function at time 't' and % a small increment, 'dt', is the probability that an event occurs% sometime between 't' and 't+dt'.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 3.  Gamma-distributed interevent intervals.% In this section we're going to simulate a renewal process by% a different method.  Imagine that you have a neuron that integrates% excitatory post-synaptic potentials that arrive as a Poisson% process.  Further, let's say that once "n" of these EPSPs have% been integrated the neuron fires a spike and resets its membrane% voltage to its resting state.  This is an example of a renewal% process: interevent intervals are independent and come from a % single underlying distribution.  We haven't explictly stated% the parametric form of the interevent interval distribution, but% that's OK.  It'll turn out that the interevent interval distribution% has a nice analytical solution.  Lucky us!% To start, let's say that the neuron has to receive exactly 3% EPSPs to fire a spike.  We'll further assume that EPSPs arrive% as a Poisson process and we'll adjust the rate of these Poisson% processes so that our model neuron will fire, on average, 50 spikes% per second (if should be clear that if we don't adjust the rate% of the EPSP processes our model neuron will fire progressively% more slowly as we require that it integrate more and more EPSPs to % fire a spike).nstepstothresh = 5;spike_rate = 50;mean_int = 1000/spike_rate;intervals = exprnd(meanint/nstepstothresh,nstepstothresh,1000);intervals = sum(intervals,1);% Let's take a look at the histogram of interspike intervals and the% spike rasters from the first 500 ms of the simulation.plot2ras(2,intervals,nstepstothresh);% Notice that we don't get many very short interspike intervals.% This is because our neuron has to wait for the 5th EPSP to fire% a spike.  In principle, 5 EPSPs could arrive immediately after our% neuron fires a spike, in which case it would fire a second spike,% but this is exceedingly unlikely to happen.% For comparison, let's simulate spikes from a neuron that only needs% to integrate 1 EPSP before it fires a spike.  It should come as % no surprise that this neuron fires as a Poisson process -- it waits % for an EPSP (which arrives as a Poisson process by assumption) and% then fires a spike.nstepstothresh = 1intervals = exprnd(meanint/nstepstothresh,nstepstothresh,1000);intervals = sum(intervals,1);plot2ras(1,intervals,nstepstothresh);% By now the fact that the interspike interval histogram looks % exponential should come as no surprise.% Notice that we occasionally get very short (and very % long) interspike intervals.  It's a lot more likely for pairs% of Poisson events to occur nearly simultaneously than for a% packet of five Poisson events to occur nearly simultaneously.% Likewise, pairs of Poisson events may have a substantial% delay between them, whereas the time between the fifth% Poisson event and the 10th (assuming that the events are coming% 5x as quickly) is unlikely to be as long.% Finally, let's look at a neuron that needs to integrate 20 Poisson% EPSPs before firing a spike.nstepstothresh = 20intervals = exprnd(meanint/nstepstothresh,nstepstothresh,1000);intervals = sum(intervals,1);plot2ras(3,intervals,nstepstothresh);% Now compare vertically across the panels.  As we require that our neuron% integrate progressively greater numbers of EPSPs, the interspike interval% histogram gets progressively tighter (less variable) and the spike trains% become progressively more regular looking.% The CV of the interspike interval distribution decreases as we% increase the number of EPSPs to count.  As was shown above, this has % the consequence of decreasing the mean to variance ratio of the spike% counts.  We fixed the mean spike count in this simulation, so we see% a decrease in the spike count variance.% This can be a problem for the simple integrate and fire model of % neuronal firing.  It predicts much less variable responses than one% actually observes in in vivo cortical recordings.  There are ways% of getting around this, however.  For instance, the fact that real% neurons receive both excitation *and inhibition* changes the situation% dramatically.% (A final aside: this is a bit tangential to the focus of this tutorial,% but it's interesting to note that the interspike interval distributions% for the integration model has a gamma density.  I'll superimpose some% gamma density functions on the interspike interval histograms to try to % convince you of this.)subplot(3,2,1); hold on;plot([0:60],gampdf([0:60],1,meanint),'g-');subplot(3,2,3); hold on;plot([0:60],gampdf([0:60],5,meanint/5),'g-');subplot(3,2,5); hold on;plot([0:60],gampdf([0:60],20,meanint/20),'g-');% (This is a nice result because it means that we could have simulated% this process simply by drawing interspike intervals from a gamma % distribution.)