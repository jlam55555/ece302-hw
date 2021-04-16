Consider a radar detection system, where we want to make a determination if a target is present or absent. If the target is present, Y = A + X, where A is a known constant. If the target is not present, Y = X. X is a zero mean Gaussian with variance σ2. The a priori probability that the target is not present is .8. The signal to noise ratio is the ratio of A to σ2, you must choose values for this exercise that give insightful results.

a) Derive and implement in MATLAB the MAP rule for detecting the target. Run 1000 iterations of your detector; compare the probability of error with the theoretical probability of error.

b) Implement a simulation that plots the receiver operating curve for this detector. Plot the receiver operating curve for several signal to noise ratios.

c) Assume that missing the target is 10 times worse than falsely detecting the target. What is the decision rule that minimizes the conditional risk? Mark this point on your receiver operating curve for at least one SNR value.

d) Using the cost structure in part c), Select one SNR value and plot the value of the expected cost for a range of a priori target present probabilities from 0 to 1.

e) Now, repeat parts a and b, but change the model such that the target present remains Y = A+X but the target not present model is now Y = A+Z where Z is a zero mean Gaussian random variable with σ2z > σ2 . Plot a few receiver operating curves for different ratios of σ2z to σ2.