referee response
==========
Major Issue I: Reorganization of Section 4
------------------------------------------
Section 4 has been reorganized in the manner suggested by the reviewer. The authors agree this structure will be more clear to the readers. 

------------------------------------------
Major Issue II: Statistical Model
------------------------------------------

We thank the referee for the very thorough analysis of this method and for the detailed suggestions. We used this input to review and revise the fitting of the ferengi data and construction of the debiasing equation for the Hubble data, as follows:

We tested a variety of models for fitting f/f_o, including linear, exponential, quadratic, and cubic, each using both a constant zeta parameter and a linear with surface brightness dependence. We compared the relative quality of fits using the AIC parameter, and found no great difference in the choice of model. The details and code can be viewed in this python notebook https://github.com/willettk/gzhubble/blob/master/python/notebooks/zeta_models_compared.ipynb. Of these, we believe the best model is an exponential fit for f_features vs z with a linear zeta, based on the following:

1) While the ferengi data could be fit relatively well using any of the models tested, the exponential is ideal because it is the only form which allows f_features to be simply and continuously bound between 0 and 1, which is a requirement of the debiased vote fractions. Further, the AIC test showed the exponential fit was either better than, or had similar support to the other models tested. 


2) We choose a linear equation to represent zeta as a function of surface brightness. In terms of goodness of fit, the AIC test showed both models were of equal quality. While this and the small slope in zeta shows we were not able to detect a surface brightness dependence, we still opt to keep the more general, linear form. As similar debiasing methods are planned for future GZ projects with a greater quantity of ferengi data, in which a possible surface brightness dependence could be better detected, keeping the general form for this project will allow for a closer comparison of the two iterations of the method. We added an emphasis on this reasoning in the paper so as not to lead the readers to the conclusion that we are claiming a strong surface brightness depedence detection in this run. 


Additionally, we chose to modify equation 2 to match the normalization of equation 4 (as suggested in minor comment 4).  

-----------------------------------
Minor Comment 1
-----------------------------------


-----------------------------------
Minor Comment 2
-----------------------------------

The FERENGI sample is biased to lower surface brightness than the real data, simply because the local Universe does not include galaxies which are as high surface brightness as are found at high redshift. We get around this partially by applying the evolutionary correction described in Section 2.3.2, but this still results in a sample which only spans, not reproduces the surface brightness-redshift distribution of the real data. This is fine as we are not correcting the population as a whole, but only galaxies at a given surface brightness and redshift.

-----------------------------------
Minor Comment 3
-----------------------------------


We have implemented changes in the text that we believe removes the confusion.

-----------------------------------
Minor Comment 4
-----------------------------------

We originally left equation 2 un-bounded because the data being fit was inherently bounded between 0 and 1, while the correction equation 4 would produce vote fractions greater than 1 if left in the same form as 2. We recognize now that the difference in normalization, while being only a small effect on the fitting of the data, is inconsistent and confusing to the reader. Therfore we changed equation 2 to match equation 4 (see response to Major Comment II). 


----------------------------------
Minor Comment 5
----------------------------------

This mistake is edited in the text.

---------------------------------
Minor Comment 6
----------------------------------

***response pending final fit ***

----------------------------------
Minor Comment 7
----------------------------------

The regions can be perfectly separated if we keep the discrete binning as in Figure 9 - however we assumed that a true/physical separation would be a smooth function of z-mu-f. Implementing a convex hull in this space created a small amount of contamination around the edges, which was minimized by adjusting the hull as described. 

***does this sound sufficient?*** 

----------------------------------
Minor Comment 8
----------------------------------

We chose to evaluate the distribution of f_features in 5 (not 4, this was a typo) bins to match the common ranges of vote fractions used in GZ morphological studies: 0-0.2 indicating a weak/no detection of the feature in question, 0.2-0.4, 0.4-0.6, and 0.6-0.8 representing different ranges of intermediate detection, and 0.8-1.0 representing a very strong detection. This explanation has been added to the text. 

----------------------------------
Minor Comment 9
----------------------------------
80% was a mostly arbitrary choice to represent the distribution of f_features in a given bin; we discussed and agree with the referee that a more typical statistical range would be more appropriate. We therefore chose to represent the range as the upper and lower 1 sigma limits, or the inner 68% of the data. 

----------------------------------
Minor Comment 10
----------------------------------
We agree that the normalized chi squared is not the ideal metric for evaluating the goodness of fit here. Rather than comparing the fits of the polynomials for Task 1 vs the rest, we chose to compare the available data per bin. For Task 1, an average of 59 points were fit in each z-mu bin (Figure 9). For the other tasks, an average of 9 points remained in each bin, due to these questions often not being asked of at least 10 users per galaxy, as is required by default in Task 1. We believe this simpler evaluation will make it more evident to the reader that the ferengi data was not well-fit for questions beyond Task 01. Last, we felt that stating these numbers was sufficient and opted to remove the example visual from the appendix. 

----------------------------------
Minor Comment 11
----------------------------------
Typo has been fixed in the text. 


