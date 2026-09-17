# CATER-MR Lean v0.9 mathematical specification

## 1. Scope

CATER-MR v0.9 is intentionally conservative. It does not introduce a new Mendelian-randomization estimator. Its core contribution is **outcome-independent qualification of trans-eQTL instruments using a cell-type-specific one-hop TF→target GRN plus full cis-eQTL evidence for the proposed parent TF**.

The primary causal estimate remains cis MR whenever a usable cis instrument set exists. Qualified trans instruments provide triangulation and optional precision augmentation.

The method is designed for:

\[
\text{full cis summary} + \text{significance-censored trans summary} + \text{directed TF→target GRN}.
\]

## 2. Notation

For target gene \(X\):

- \(Pa(X)\): direct one-hop TF parents of \(X\) in the cell-type GRN;
- \(\hat\beta_{GX}\): SNP→target expression effect;
- \(\hat\beta_{GT}\): the same SNP→parent-TF cis-eQTL effect;
- \(\hat\beta_{GY}\): SNP→outcome effect;
- \(R\): signed, allele-aligned LD correlation matrix;
- \(D_X\), \(D_Y\): diagonal matrices of exposure and outcome standard errors.

For retained instruments:

\[
\gamma=(\hat\beta_{G_1X},\ldots,\hat\beta_{G_KX})^T,
\qquad
\Gamma=(\hat\beta_{G_1Y},\ldots,\hat\beta_{G_KY})^T.
\]

## 3. Cis instruments

The target cis instrument set is

\[
\mathcal C_X
=
\{G:\;G\in cis(X),\;P_{GX}<\tau_X\},
\]

with default

\[
\tau_X=5\times10^{-8}.
\]

## 4. Raw trans candidates

The GRN is used only as a one-hop topology prior. A reported trans association is a raw candidate when

\[
\mathcal T_X^{raw}
=
\{G:\;P_{GX}^{trans}<\tau_X,\;G\in cis(T),\;T\in Pa(X)\}.
\]

Physical proximity to a parent TF is **candidate generation only**. It is not sufficient evidence that the SNP perturbs that TF.

## 5. Parent-TF cis anchoring

For each \(G\in\mathcal T_X^{raw}\), CATER-MR queries the full cis summary statistics of every direct parent TF whose cis interval contains \(G\).

A parent TF qualifies when

\[
P_{GT}<\tau_T,
\]

with default

\[
\tau_T=\tau_X=5\times10^{-8}.
\]

Let

\[
Q_X(G)=\{T\in Pa(X):G\in cis(T),\;P_{GT}<\tau_T\}.
\]

The mechanism-supported trans set requires a unique qualifying parent:

\[
\mathcal T_X^{extended}
=
\{G\in\mathcal T_X^{raw}:|Q_X(G)|=1\}.
\]

If \(|Q_X(G)|>1\), the SNP is marked `TRANS_AMBIGUOUS_PARENT` and excluded from the default trans set because the proposed source regulator is not identifiable from the available evidence.

The supported mechanism is therefore

\[
G\xrightarrow{cis\ eQTL}T
\xrightarrow{GRN}X
\]

with an independently observed

\[
G\xrightarrow{trans\ eQTL}X.
\]

This increases mechanistic plausibility but does **not** prove mediation or the MR exclusion restriction.

## 6. Observed trans-specificity filter

For significance-censored trans catalogs, define

\[
H_G
=
\left|\{Z:(G,Z)\text{ is a reported significant trans-eQTL pair}\}\right|.
\]

The default core trans set is

\[
\mathcal T_X^{core}
=
\{G\in\mathcal T_X^{extended}:H_G\le h_{max}\},
\]

with default

\[
h_{max}=1.
\]

This is only an **observed significant-target burden**. If another SNP-gene association is absent from a significance-censored catalog, CATER-MR treats it as unavailable/censored information, not as \(\beta=0\).

`trans_set="extended"` removes the \(H_G\) filter but retains the TF-cis anchor and unique-parent requirements.

## 7. LD selection

Candidates are ordered by target-association P value. A candidate \(G_j\) is retained only if for every already-retained \(G_k\) within the configured genomic window,

\[
r_{jk}^2<\tau_{LD}.
\]

Default:

\[
\tau_{LD}=0.01.
\]

Residual signed LD among retained variants is preserved for MR estimation.

## 8. Generalized IVW estimator

For any retained instrument set,

\[
\Omega_Y=D_YRD_Y.
\]

For \(K>1\):

\[
\hat\theta
=
\frac{\gamma^T\Omega_Y^{-1}\Gamma}
     {\gamma^T\Omega_Y^{-1}\gamma},
\]

\[
SE(\hat\theta)
=
\left(\gamma^T\Omega_Y^{-1}\gamma\right)^{-1/2}.
\]

Residual heterogeneity:

\[
Q
=
(\Gamma-\hat\theta\gamma)^T
\Omega_Y^{-1}
(\Gamma-\hat\theta\gamma),
\qquad df=K-1.
\]

For one IV:

\[
\hat\theta=\frac{\Gamma}{\gamma},
\qquad
SE(\hat\theta)=\left|\frac{SE_Y}{\gamma}\right|.
\]

## 9. Core models

CATER-MR reports three substantive fits:

### Cis primary

\[
\hat\theta_{cis}=GIVW(\mathcal C_X).
\]

If available, this is always the primary model.

### Qualified trans

\[
\hat\theta_{trans}=GIVW(\mathcal T_X^{core})
\]

for the default analysis, or \(\mathcal T_X^{extended}\) when explicitly requested.

This is a triangulation/sensitivity estimate.

### Augmented CATER

\[
\hat\theta_{aug}
=
GIVW(\mathcal C_X\cup\mathcal T_X^{core}).
\]

The implementation retains a historical `combined` alias for backward compatibility, but documentation uses `augmented`.

The augmented fit is **not automatically promoted to primary**, even if it is more precise.

## 10. Why precision can improve without validity improving

If the structural relation is

\[
\Gamma=\theta\gamma+\alpha+\varepsilon,
\]

where \(\alpha\) contains SNP→outcome effects not mediated through target \(X\), then

\[
E(\hat\theta)
=
\theta+
\frac{\gamma^T\Omega_Y^{-1}\alpha}
     {\gamma^T\Omega_Y^{-1}\gamma}.
\]

Adding valid trans IVs increases information, but adding trans IVs with structured pleiotropy can increase bias. Therefore lower standard error is not interpreted as proof of a better causal estimate.

## 11. Instrument strength

Exposure covariance:

\[
\Sigma_X=D_XRD_X.
\]

Joint exposure information:

\[
I_X=\gamma^T\Sigma_X^{-1}\gamma.
\]

Effective strength:

\[
F_{eff}=\frac{I_X}{K}.
\]

v0.9 reports separate values for cis, trans, augmented, and primary models. The backward-compatible generic `effective_F` field maps to the actual primary model.

## 12. Precision gain

Define

\[
J=\gamma^T\Omega_Y^{-1}\gamma.
\]

The trans contribution to augmented MR precision is

\[
f_{trans,prec}
=
\frac{J_{aug}-J_{cis}}{J_{aug}}.
\]

The directly interpretable standard-error reduction is

\[
R_{SE}
=
1-\frac{SE_{aug}}{SE_{cis}}.
\]

These metrics quantify precision gain only; they do not measure reduction of pleiotropic bias.

## 13. Cis-trans heterogeneity

Let \(w_c\) and \(w_t\) be normalized GIVW weight vectors for cis and trans subsets and

\[
\Omega_{ct}=D_{Y,c}R_{ct}D_{Y,t}.
\]

Then

\[
Cov(\hat\theta_c,\hat\theta_t)
=w_c^T\Omega_{ct}w_t,
\]

\[
Z_{c-t}
=
\frac{\hat\theta_c-\hat\theta_t}
{\sqrt{SE_c^2+SE_t^2-2Cov(\hat\theta_c,\hat\theta_t)}}.
\]

This is a consistency diagnostic, not an outcome-based IV-selection rule.

## 14. TF-locus influence diagnostics

For a parent-TF locus \(T\):

\[
\Delta_T
=
|\hat\theta_{all}-\hat\theta_{-T}|.
\]

Exposure-information weight:

\[
w_{T,exp}
=
\frac{I_{all}-I_{-T}}{I_{all}}.
\]

MR-precision weight:

\[
w_{T,prec}
=
\frac{J_{all}-J_{-T}}{J_{all}}.
\]

These diagnose whether one parent locus dominates the augmented result.

## 15. Advanced sensitivity analyses

Sibling co-perturbation screening and correlated-IV MVMR remain implemented but are disabled by default.

A standard MVMR fit requires a complete selected-SNP × exposure effect matrix \(B\):

\[
\hat\theta_{MV}
=
(B^T\Omega_Y^{-1}B)^{-1}B^T\Omega_Y^{-1}\Gamma_Y.
\]

If any required cross-exposure SNP effect or standard error is unavailable, the implementation fails closed with `MVMR_CROSS_EXPOSURE_EFFECTS_UNAVAILABLE`. Missing associations from a significant-only trans catalog are never filled with zero.

The correlated conditional-F implementation remains explicitly experimental and cannot make a network model primary in Lean v0.9.

## 16. Claims the method does not make

The following implications are explicitly rejected:

\[
G\in locus(T) \not\Rightarrow G\rightarrow T,
\]

\[
G\rightarrow T\;\&\;T\rightarrow X\;\&\;G\rightarrow X
\not\Rightarrow
\text{exclusion restriction holds},
\]

\[
\text{not reported trans association}\not\Rightarrow\beta=0,
\]

\[
SE_{aug}<SE_{cis}\not\Rightarrow\hat\theta_{aug}\text{ is less biased}.
\]

The GRN and TF-cis anchor are therefore interpreted as **mechanism filters that reduce, but cannot eliminate, trans-IV validity risk**.

## References

1. Aguet F, Brown AA, Castel SE, et al. Genetic effects on gene expression across human tissues. *Nature*. 2017;550:204-213. https://doi.org/10.1038/nature24277
2. Zheng J, Haberland V, Baird D, et al. Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. *Nature Genetics*. 2020;52:1122-1131. https://doi.org/10.1038/s41588-020-0682-6
3. Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC. Mendelian randomization with fine-mapped genetic data: Choosing from large numbers of correlated instrumental variables. *Genetic Epidemiology*. 2017;41:714-725. https://doi.org/10.1002/gepi.22077
4. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *International Journal of Epidemiology*. 2019;48:713-727. https://doi.org/10.1093/ije/dyy262
