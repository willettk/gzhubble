referee response
==========
Major Issue I: Reorganization of Section 4
------------------------------------------
Section 4 has been reorganized in the manner suggested by the reviewer. The authors agree this structure will be more clear to the readers. 

------------------------------------------
Major Issue II: Statistical Model
------------------------------------------

We thank the referee for the very thorough analysis of this method and for the detailed suggestions. We used this input to review and revise the fitting of the ferengi data and construction of the debiasing equation for the Hubble data, as follows:

We tested a variety of models for fitting f/f_o, including linear, exponential, polynomial, and power, each testing a constant zeta parameter and a linear with surface brightness dependence. We compared the relative quality of fits using the AIC parameter, and found no great difference in the choice of model (***should quantify this if we can***), the details and code can be viewed in this python notebook (***link eventually***). Of these, we believe the best model is an exponential fit for f_features vs z with a linear/constant (***need to choose***) zeta, based on the following:

1) While the ferengi data could be fit relatively well using any of the models tested, the exponential is ideal because it is the only form which allows f_features to be simply and continuously bound between 0 and 1, which is a requirement of the debiased vote fractions. 

***Zeta: need to choose which.***

2a) We choose a linear equation to represent zeta as a function of surface brightness because (***I can't think of a great reason other than it doesn't hurt since the slope is ~0 anyway, but that feels like such a weak argument...***) 

2b) We choose to remove the surface brightness dependence from zeta and use a constant parameter. Between the exponential fit with a constant and linear model for zeta, the constant model was shown to be a slightly better fit via the AIC parameter. We note that this does not change the effect of the debiasing equation strongly, as the original slope was approximately zero. This change is mostly implemented to remove the possible implication that there is any strong surface brightness dependence in the paper. 

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

(*** will continue to be taken care of as section 4 continues to be finalized (pending choice of zeta in Major Comment II), but probable response will be:***)

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

***response pending final fit / choice of zeta***

----------------------------------
Minor Comment 7
----------------------------------

The regions can be perfectly separated if we keep the discrete binning as in Figure 9 - however we assumed that a true/physical separation would be a smooth function of z-mu-f. Implementing a convex hull in this space created a small amount of contamination around the edges, which was minimized by adjusting the hull as described. 

***does this sound efficient?*** 

----------------------------------
Minor Comment 8
----------------------------------



