referee response
==========
Major Issue I: Reorganization of Section 4
------------------------------------------
Section 4 has been re-organized in the manner suggested by the reviewer. The authors agree this structure will be more clear to the readers. 

------------------------------------------
Major Issue II: Statistical Model
------------------------------------------

We thank the referee for the very thorough analysis of this method and for the detailed suggestions. We used this input to review and revise the fitting of the artificially-redshifted images (from FERENGI) and construction of the debiasing equation for the Hubble data, as follows:

We tested a variety of models for fitting f/f_0, including linear, quadratic, cubic, and exponential dependences on redshift; for each of these functional forms, we also tested combinations using either a constant zeta parameter or a zeta function that linearly depends on surface brightness. We compared the relative quality of fits using the Akaike information criteria (AIC); for each case, we found no significant differences in the results based on the choice of model. Full details and code for these tests can be viewed in this Python notebook: https://github.com/willettk/gzhubble/blob/master/python/notebooks/zeta_models_compared.ipynb. 

Of all the models tested, we believe the best model is an exponential fit for f_features vs redshift using a linear dependence on zeta. This decision is based on the following points:

1) While the FERENGI data can be fit relatively well using any of the above models, the exponential is ideal because it is the only form which allows f_features to be simply and continuously bound between 0 and 1. This is a key requirement for the debiased votes so that they can be interpreted as normalized fractions. Furthermore, the results of the AIC test show that the exponential fit is either better than, or had similar support to the other models tested. 

2) We choose a linear equation to represent zeta as a function of surface brightness. In terms of goodness-of-fit, the AIC test showed both models were of equal quality. While this and the small slope in zeta shows that there is no significant surface brightness dependence in the GZH data, we still opt to keep the more general, linear form. As similar debiasing methods are planned for future GZ projects with a greater quantity of calibration images (and for which there may exist stronger dependencies on surface brightness), keeping the general form for this project will allow for a closer comparison of the two iterations of the method. We have added an emphasis on this reasoning in the manuscript so as not to lead the readers to the conclusion that we are claiming a strong surface brightness depedence detection in this run. 

Finally, we have modified Equation 2 to match the normalization of Equation 4 (as suggested in minor comment 4).  

-----------------------------------
Minor Comment 1
-----------------------------------

***1) Section 3.2: what checking did the authors do to indicate that giving all anonymous classifications a weight of 1 is valid?***

From a strict sense, it is impossible to tell a priori if giving the anonymous classifiers a weight of 1 is valid (since the grouping within the anonymous data is unknown). However, there are large amounts of data from similar projects that support this approach. For the Galaxy Zoo: CANDELS project (Simmons et al. 2016, submitted), which is nearly identical in design to Galaxy Zoo: Hubble and shares a large fraction of the same user base, the distribution of the classifier consistency for both the logged-in and non-logged-in populations peak at the same value (kappa = 0.7). The anonymous users have slightly broader wings, with slightly more addresses at the lowest and highest levels of consistency, but the overall behavior is similar. The roughly 50-50 split between logged-in users and non-logged-in IP addresses is typical of the ratio for Zooniverse crowdsourced projects as a group.

ATTACH: Figure from https://cloud.githubusercontent.com/assets/2357126/17576650/cfc59fbe-5f29-11e6-99b7-c39a43c013b6.png

-----------------------------------
Minor Comment 2
-----------------------------------

***2) Figure 6: the authors should comment (unless I missed the comment) on why the bivariate distribution of the FERENGI galaxies in (mu,z) does not match the bivariate distribution for the main sample.***

The FERENGI sample is biased in including examples of galaxies that have lower surface brightnesses than any images found in the real Hubble data, simply because the local Universe has no examples of galaxies at the high surface brightness end. This is a combination of both a finite-volume sample at low redshift in addition to some evolution in the stellar populations and size-magnitude relation over cosmic time.  We deal with this issue partially by applying the evolutionary correction described in Section 2.3.2, but this still results in a sample spans (but does not reproduce) the surface brightness-redshift distribution of the real data. The mismatch does not affect the overall calibration accuracy as we are not correcting the population as a whole, but only galaxies at particular bins of surface brightness and redshift.

We have added a paragraph in Section 4 for clarification.

-----------------------------------
Minor Comment 3
-----------------------------------

***3) In parts of Section 4.1, there is a confusing mix of f_mu,z and f_features, which appear to represent the same thing. Some clean-up (or additional explanation as to why they are different quantities) would be useful here.***

The main difference between Equations 2 and 4 is that the former is applied to the FERENGI calibration images and used to fit the normalization and zeta parameters. The latter uses the parameters to apply corrections to individual galaxies from the true data, and is the form implemented throughout the text. 

We have added text to explicitly mark the difference between the two equations in what is now Section 4.2. 

-----------------------------------
Minor Comment 4
-----------------------------------

***4) I probably should have mentioned this above, but why are the authors not working with equation (4) as the basis of their regression model?***

Equation 2 was originally left unbounded because the raw data being fit was inherently bounded between vote fractions of 0 and 1. However, the original verseion of Equation 4 could produce vote fractions greater than 1 if left in the same form as Equation 2. We have realized that the difference in normalization, while having only a small effect on the fitting of the data, is inconsistent and may be confusing to a reader. Therefore, we have changed Equation 2 to match Equation 4 (see response to Major Comment II, as well as Minor Comment 3). All data in the tables and catalog have been updated to account for these small changes.

----------------------------------
Minor Comment 5
----------------------------------

***5) Figure 8: "-10 < log(xi) < 10" should be "-10 < xi < 10"***

This mistake is edited in the text.

---------------------------------
Minor Comment 6
----------------------------------

***6) How many data points are there in Figure 8? Obviously more than 37, but fewer than 37 x 8.***

***response pending final fit ***

----------------------------------
Minor Comment 7
----------------------------------

***7) Section 4.2: "The boundaries [of the convex hull] are then adjusted until the contamination from both groups is minimized." If I look at Figure 9, it would appear that the correctable and lower-limit regions can be perfectly separated (using a tree, for instance). Am I wrong about this?***

The regions in Figure 9 can be perfectly separated in f_features-redshift space if we keep the same discrete binning shown in the figure. However, our analysis makes the assumption that a true/physical separation should be a smooth function of z-mu-f. Implementing a convex hull in this space created a small amount of contamination around the edges, which was minimized by adjusting the hull as described. 

***does this sound sufficient?*** 

----------------------------------
Minor Comment 8
----------------------------------

***8) Section 4.2: "In each unclassifiable (shaded blue) z,mu bin, the spread of intrinsic values for f_features,z=0.3 for four ranges of observed f_features is computed." Where did the "four" come from? More explanation needed.***

We chose to evaluate the distribution of f_features in 5 (not 4; a typo in the earlier manuscript has been corrected) bins to match the common ranges of vote fractions used in previous GZ morphological studies: 0.0-0.2 is typically used to indicate weak/no detection of the feature in question, 0.2-0.4, 0.4-0.6, and 0.6-0.8 represent ranges of increasing likelihood of detection, and 0.8-1.0 represents a very strong detection. This explanation and references to other papers using a similar system has been added to the text. 

----------------------------------
Minor Comment 9
----------------------------------

***9) Why 80% of the data (orange bars, Figure 10)? This is an ad hoc choice that doesn't correspond to any typical statistical practice (except perhaps a 10% trimmed mean, but the authors are not computing means here). And how is this information applied in Figure 11, middle panel? This is not clear.***

The threshold of 80% was a mostly arbitrary choice to represent the distribution of f_features in a given bin; the authors have discussed it and agree with the referee that a more typical statistical range would be more appropriate. We therefore have changed the range to upper and lower 1 sigma limits, or the inner 68% of the data. 

----------------------------------
Minor Comment 10
----------------------------------

***10) Section 4.3: "The goodness-of-fit was evaluated in each bin, for each Task, using a normalized chi^2 metric." The authors need to provide the formula for this; the numbers don't make sense, if one assumes "normalized" is the same
as "reduced chi^2" (where the typical value for a good fit is approximately1).***

We agree that the normalized chi-squared is not the ideal metric for evaluating the goodness-of-fit here. Rather than comparing the fits of the polynomials for Task 1 vs the rest of the tasks, we compare the available data per bin. For Task 1, an average of 59 points were fit in each z-mu bin (Figure 9). For the other tasks, an average of 9 points remained in each bin, due to these questions often not being asked of at least 10 users per galaxy (as required by default in Task 1). This simpler evaluation should make it more evident to the reader that the FERENGI data was not well-fit for questions beyond Task 01. 

The authors feel that quoting the numbers is sufficient for the manuscript; we have removed the example visual from the Appendix for streamlining. 

----------------------------------
Minor Comment 11
----------------------------------

***11) Section 7.3: "external pulblications" -> "external publications"***

Typo has been fixed in the text. 


