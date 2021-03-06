We thank the referee for the very thorough analysis and for the detailed
suggestions that helped to improve the analysis and the paper overall. In the
following, we provide the original comment by the referee, followed by our
response. 

The major changes in the manuscript have been marked in red boldface. We note
that the reorganization of Section 4 (Major Issue I) is not marked. There has
also been a mild amount of reordering of tables and figures as a result. We
have added a single figure (Figure 11) that shows the fraction of correctable
galaxies distributed over the sample, and a short amount of explanatory text in
Sections 4.2 and 8. 

========== Major Issue I: Reorganization of Section 4
------------------------------------------

***I spent an inordinate amount of time trying to piece together the analysis
of Section 4 before discovering that important information was given later (in
Section 4.2) rather than earlier (in Section 4.1).***

***I believe that Section 4.2 should come first: present the FERENGI sample,
and discuss the scheme used to select those data that are to be analyzed to
determine the debiasing correction. Then the details of the statistical model
for debiasing (currently in Section 4.1) should come next.***

***Section 4.3 is fine where it is placed. (However, see minor comment 10
below, which actually is not quite as minor as the other comments.)***

Section 4 has been re-organized in the manner suggested by the reviewer. The
authors agree this structure will be more clear to the readers. 

========== Major Issue II: Statistical Model
------------------------------------------

***Due to the length of Major Issue II, it is placed at the end of this
response.*** 

We used the referee's input to review and revise the fitting of the
artificially-redshifted images (from FERENGI) and construction of the debiasing
equation for the Hubble data, as follows:

We tested a variety of models for fitting f/f_0, including linear, quadratic,
cubic, and exponential dependences on redshift; for each of these functional
forms, we also tested combinations using either a constant zeta parameter or a
zeta function that linearly depends on surface brightness. We compared the
relative quality of fits using the Akaike information criteria (AIC); for each
case, we found no significant differences in the results based on the choice of
model. Full details and code for these tests can be viewed in this Python
notebook:

https://github.com/willettk/gzhubble/blob/master/python/notebooks/zeta_models_compared.ipynb 

Of all the models tested, we believe the best model is an exponential fit for
f_features vs redshift using a linear dependence on zeta. This decision is
based on the following points:

1) While the FERENGI data can be fit relatively well using any of the above
models, the exponential is ideal because it is the only form whose limiting
behaviour can be set to approach 0 and 1. This is a key requirement for the
debiased votes so that they can be interpreted as normalized fractions.
Furthermore, the results of the AIC test show that the exponential fit is
either better than, or had similar support to the other models tested. 

2) We choose a linear equation to represent zeta as a function of surface
brightness. In terms of goodness-of-fit, the AIC test showed both models were
of equal quality. While this and the small slope in zeta shows that there is no
significant surface brightness dependence in the GZH data, we still opt to keep
the more general, linear form. As similar debiasing methods are planned for
future GZ projects with a greater quantity of calibration images (and for which
there may exist stronger dependencies on surface brightness), keeping the
general form for this project will allow for a closer comparison of the two
iterations of the method. We have added an emphasis on this reasoning in the
manuscript so as not to lead the readers to the conclusion that we are claiming
a strong surface brightness depedence detection in this run. 

Finally, we have modified Equation 2 to match the normalization of Equation 4
(as suggested in minor comment 4).  

----------------------------------- Minor Comment 1
-----------------------------------

***1) Section 3.2: what checking did the authors do to indicate that giving all
anonymous classifications a weight of 1 is valid?***

From a strict sense, it is impossible to tell a priori if giving the anonymous
classifiers a weight of 1 is valid (since the grouping within the anonymous
data is unknown). However, there are large amounts of data from similar
projects that support this approach. For the Galaxy Zoo: CANDELS project
(Simmons et al. 2016, submitted), which is nearly identical in design to Galaxy
Zoo: Hubble and shares a large fraction of the same user base, the distribution
of the classifier consistency for both the logged-in and non-logged-in
populations peak at the same value (kappa = 0.7). The anonymous users have
slightly broader wings, with slightly more IP addresses at the lowest and
highest levels of consistency, but the overall behavior is similar. 

The roughly 50-50 split between logged-in users and non-logged-in IP addresses
is typical of the ratio for Zooniverse crowdsourced projects as a group.
Similar results have been published for the Planet Hunters project (Schwamb et
al. 2012) and have been seen in unpublished data for Pulsar Hunters
(https://www.zooniverse.org/projects/zooniverse/pulsar-hunters). 

ATTACH: Figure from
https://cloud.githubusercontent.com/assets/2357126/17576650/cfc59fbe-5f29-11e6-99b7-c39a43c013b6.png

----------------------------------- Minor Comment 2
-----------------------------------

***2) Figure 6: the authors should comment (unless I missed the comment) on why
the bivariate distribution of the FERENGI galaxies in (mu,z) does not match the
bivariate distribution for the main sample.***

The FERENGI sample is biased in including examples of galaxies that have lower
surface brightnesses than any images found in the real Hubble data, simply
because the local Universe has no examples of galaxies at the high surface
brightness end. This is a combination of both a finite-volume sample at low
redshift in addition to some evolution in the stellar populations and
size-magnitude relation over cosmic time.  We deal with this issue partially by
applying the evolutionary correction described in Section 2.3.2, but this still
results in a sample that spans (but does not reproduce) the surface
brightness-redshift distribution of the real data. The mismatch does not affect
the overall calibration accuracy as we are not correcting the population as a
whole, but only galaxies at particular bins of surface brightness and redshift.

We have added a paragraph in Section 4 for clarification.

----------------------------------- Minor Comment 3
-----------------------------------

***3) In parts of Section 4.1, there is a confusing mix of f_mu,z and
f_features, which appear to represent the same thing. Some clean-up (or
additional explanation as to why they are different quantities) would be useful
here.***

The main difference between Equations 2 and 4 is that the former is applied to
the FERENGI calibration images and used to fit the normalization and zeta
parameters. The latter uses the parameters to apply corrections to individual
galaxies from the true data, and is the form implemented throughout the text. 

We have added text to explicitly mark the difference between the two equations
in what is now Section 4.2. 

----------------------------------- Minor Comment 4
-----------------------------------

***4) I probably should have mentioned this above, but why are the authors not
working with equation (4) as the basis of their regression model?***

Equation 2 was originally left unbounded because the raw data being fit was
inherently bounded between vote fractions of 0 and 1. However, the original
version of Equation 4 could produce vote fractions greater than 1 for certain
combinations of z and f_0 if left in the same form as Equation 2. We have
realized that the difference in normalization, while having only a small effect
on the fitting of the data, is inconsistent and may be confusing to a reader.
Therefore, we have changed Equation 2 to match Equation 4 (see response to
Major Comment II, as well as Minor Comment 3). All data in the tables and
catalog have been updated to account for these small changes.

The attached figure displays the effect of this change. Shown is the
distribution of the change in debiased f_features between this draft (new
normalization) and the previous draft (old normalization). For the correctable
sample of galaxies, the new debiased f_features vote fractions are on average
0.02 greater than the previous draft, and at maximum 0.06 greater. These
changes are too small to significantly affect any final morphologies. 

ATTACH:
https://github.com/willettk/gzhubble/blob/master/python/notebooks/compare_new_debiased_fraction.png

---------------------------------- Minor Comment 5
----------------------------------

***5) Figure 8: "-10 < log(xi) < 10" should be "-10 < xi < 10"***

Correct, thank you. This typo has been edited in the caption (now Figure 10).

--------------------------------- Minor Comment 6
----------------------------------

***6) How many data points are there in Figure 8? Obviously more than 37, but
fewer than 37 x 8.***

The final version of this figure (now Figure 10; changed slightly based on the
new form of Equation 4) has 28 unique galaxies, each with multiple images that
have been artificially redshifted in between three to eight bins. 

---------------------------------- Minor Comment 7
----------------------------------

***7) Section 4.2: "The boundaries [of the convex hull] are then adjusted until
the contamination from both groups is minimized." If I look at Figure 9, it
would appear that the correctable and lower-limit regions can be perfectly
separated (using a tree, for instance). Am I wrong about this?***

The regions in Figure 9 (now Figure 7) can be nearly perfectly separated in
f_features-redshift space, but there are a few regions near the edges in which
there is overlap between the two populations. This is affected slightly in the
figure by the choice of bin size and location. However, our analysis makes the
assumption that a true/physical separation should be a smooth function of
z-mu-f. Implementing a convex hull in this space created a very small volume
near the edge of the inner convex hull that had overlapping positive and
negative data points; we optimize the overall accuracy of the hull by adjusting
the outer boundaries as described in what is now Section 4.1. 

Use of a different method such as a tree could potentially prove to be a
superior separator. However, the new normalization technique used in this draft
results in very clean separation based on the convex hull, with only a single
adjustment; see:

https://github.com/willettk/gzhubble/blob/master/python/creating_debiased_catalog/STEP_3_ferengi_smooth_function.ipynb

Further optimization would have an extremely minimal effect on the overall
results. 

---------------------------------- Minor Comment 8
----------------------------------

***8) Section 4.2: "In each unclassifiable (shaded blue) z,mu bin, the spread
of intrinsic values for f_features,z=0.3 for four ranges of observed f_features
is computed." Where did the "four" come from? More explanation needed.***

We chose to evaluate the distribution of f_features in 5 (not 4; a typo in the
earlier manuscript has been corrected) bins to match some of the common ranges
of vote fractions used in previous GZ morphological studies. The bin size is
chosen to be roughly twice the size of the typical scatter for any individual
galaxy (which has been measured using repeat image classifications in several
phases of Galaxy Zoo), and so any finer grained intervals will begin to simply
fit the noise. 

Some clarification has been added to the text at the end of Section 4.1. 

---------------------------------- Minor Comment 9
----------------------------------

***9) Why 80% of the data (orange bars, Figure 10)? This is an ad hoc choice
that doesn't correspond to any typical statistical practice (except perhaps a
10% trimmed mean, but the authors are not computing means here). And how is
this information applied in Figure 11, middle panel? This is not clear.***

The authors agree with the referee that a more typical statistical range would
be appropriate. We have changed the range to upper and lower 1 sigma limits, or
the inner 68% of the data. The figure (now Figure 8) has been modified and text
added to the end of Section 4.1.

The edges of the 1-sigma range in each bin are used to determine the lower and
upper limits for the debiased vote fractions of f_features. In the middle panel
of what is now Figure 9, we show only the lower limits.

---------------------------------- Minor Comment 10
----------------------------------

***10) Section 4.3: "The goodness-of-fit was evaluated in each bin, for each
Task, using a normalized chi^2 metric." The authors need to provide the formula
for this; the numbers don't make sense, if one assumes "normalized" is the same
    as "reduced chi^2" (where the typical value for a good fit is approximately
    1).***

We agree that the normalized chi-squared is not the ideal metric for evaluating
the goodness-of-fit here. Rather than comparing the fits of the polynomials for
Task 1 vs the rest of the tasks, we compare the available data per bin. For
Task 1, an average of 59 points were fit in each z-mu bin (now Figure 7).
Fitting parameters for any other task is significantly limited by the number of
galaxies in each z-mu bin; the higher-order tasks have an average of only 9
galaxies per bin. This is driven by the fact that the tasks do not meet the
minimum threshold of 10 individual classifications (this threshold is exceeded
for every galaxy in Task 1, however). This simpler evaluation should make it
    more evident to the reader that the FERENGI data was not well-fit for
    questions beyond Task 1.

As this is not a key product of the paper, we have removed Figures A5 and A6
from the manuscript to avoid confusion and extraneous discussion. Additional
text has been added at the end of Section 4.3.  

---------------------------------- Minor Comment 11
----------------------------------

***11) Section 7.3: "external pulblications" -> "external publications"***

Thank you, this typo has been fixed in the text. 

------------------------------- ***Major Issue II***
----------------------------------

***The authors currently apply a two-step sequential linear regression model.
The first step is to fit***

***f = f_o exp[-(z-z_o)/xi]***

***which is equivalent to the fitting a linear regression after
alog-transformation:***

***log(f) = log(f_o) + z_o/xi - z/xi = beta_0 + beta_1 z .***

***is fit. For the second linear model, the xi_1 term is not actually
statistically significant, but the authors sweep this under the rug and simply
use the point estimate as is to drive their debiasing correction. (I did
combine the steps into one linear regression model with interaction but
discovered that the model could not be fit directly. It could perhaps be fit
via iteration but I did not explore this further.)***

***The model appears to be a bit contrived, and there is no evidence presented
that the model actually generates good results [other than a vague mention of
chi^2 = 0.04 in Section 4.3 that may or may not be applicable here]). Since
there is no physical motivation for the model, I would suggest two alternatives
that are simpler and may perform just as well or better. (At a minimum I would
want to see a metric of fit defined [*] and computed for the current model and
the alternatives, if for no other reason to show that the current model is just
as good or better than the others. If it is "just as good," I'd go with a
simpler alternative.) It may seem that I am asking for too much here, but since
the whole goal of this paper is to provide a catalog with good estimates of
f_features that one can use to select galaxy samples, the authors need to do
more to show that their debiased values of f_features are actually valid!
(Estimates of uncertainty would help too; the models below would be able to
provide them.)***

***[*] Is the normalized chi^2 considered a metric of fit? It is not completely
clear. See minor comment 10 below. Perhaps the mean-squared error would work
better here; it is a standard estimator.***

***1) Linear regression, with interactions; perhaps***

***log(f/f_o) = beta_0 + beta_1 z + beta_2 mu + beta_12 z mu ,***

***although one should check other tranformations (e.g., sqrt(f/f_o), log(mu),
etc.)***

***2) Forego parametrization, and apply a regression tree (or its related
cousin random forest regression) to log(f/f_o)(mu,z). Play with the tree
settings to achieve a proper trade-off between bias and variance (i.e., to make
sure the nodes are numerous enough to capture the true dependence of log(f/f_o)
on mu and z, but not so numerous that each node only has a few galaxies, making
the estimates of log(f/f_o) unduly noisy). A good book to use to learn about
regression trees (and linear regression too) is "An Introduction to Statistical
Learning" by Gareth James et al. The authors should not worry about the
"Applications in R" part, although they might find doing the analysis in R
(assuming it is not already done there) might be simplest. Note that a PDF of
the book is available for free from James' website.***

***Another issue in Section 4.1 is the statement "The combination of all such
parameters forms a high-dimensional space, and it is not clear how to separate
this into individual effects." (1) I'm not sure why the individual effects are
important here: the goal is to derive a debiasing correction, not to understand
how the debiasing correction works. (2) There are many methods that one can use
to deal with higher-dimensional regression problems, such as PCA regression,
the Lasso, best subset selection, etc.; the book mentioned above gives more
details on these. I'm not going to demand that the authors work with the
high-dimensional data now, but they should begin to learn these statistical
methods for future work.***
