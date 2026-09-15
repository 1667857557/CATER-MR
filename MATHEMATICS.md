# CATER-MR mathematical specification

## Notation

For retained instruments \(G_1,\ldots,G_K\):

- \(\gamma=(\hat\beta_{G_1X},\ldots,\hat\beta_{G_KX})^T\): marginal SNP-exposure effects.
- \(\Gamma=(\hat\beta_{G_1Y},\ldots,\hat\beta_{G_KY})^T\): harmonized marginal SNP-outcome effects.
- \(D_X=\mathrm{diag}(SE_{G_1X},\ldots,SE_{G_KX})\).
- \(D_Y=\mathrm{diag}(SE_{G_1Y},\ldots,SE_{G_KY})\).
- \(R\): signed LD correlation matrix after allele alignment.

## Instrument eligibility

For target gene \(X\), common statistical eligibility is

\[
P_{GX}<P_{\mathrm{instrument}}.
\]

Cis candidates are

\[
\mathcal C_X^{\mathrm{cis}}
=
\{G:\;G\in\mathrm{cis}(X),\;P_{GX}<P_{\mathrm{instrument}}\}.
\]

For direct GRN parents \(Pa(X)\), trans candidates are

\[
\mathcal C_X^{\mathrm{trans}}
=
\{G:\;P_{GX}<P_{\mathrm{instrument}},\;G\in\mathrm{locus}(T),\;T\in Pa(X)\}.
\]

For significance-censored trans catalogs,

\[
P_{\mathrm{instrument}}\le P_{\mathrm{reporting}}^{\mathrm{trans}}.
\]

Any connected genomic component containing the target cis interval is assigned to the cis component.

## LD selection

Candidates are ordered by target-association \(P\)-value. A candidate \(G_j\) is retained when, for every previously retained \(G_k\) within the configured genomic window,

\[
r_{jk}^{2}<\tau_{LD}.
\]

The retained signed \(R\) is propagated to MR estimation.

## Generalized IVW

Outcome covariance:

\[
\Omega_Y=D_YRD_Y.
\]

For \(K>1\),

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

For one IV,

\[
\hat\theta=\frac{\Gamma}{\gamma},
\qquad
SE(\hat\theta)=\left|\frac{SE_Y}{\gamma}\right|.
\]

## Instrument strength

Per-IV statistic:

\[
F_j=\left(\frac{\gamma_j}{SE_{X,j}}\right)^2.
\]

Exposure covariance:

\[
\Sigma_X=D_XRD_X.
\]

Joint exposure information:

\[
I_X=\gamma^T\Sigma_X^{-1}\gamma.
\]

Reported effective strength:

\[
F_{\mathrm{eff}}=\frac{I_X}{K}.
\]

## Trans information fractions

Exposure-side trans increment:

\[
f_{\mathrm{trans,exp}}
=
\frac{I_{\mathrm{combined}}-I_{\mathrm{cis}}}
     {I_{\mathrm{combined}}}.
\]

MR-precision increment:

\[
f_{\mathrm{trans,prec}}
=
\frac{J_{\mathrm{combined}}-J_{\mathrm{cis}}}
     {J_{\mathrm{combined}}},
\qquad
J=\gamma^T\Omega_Y^{-1}\gamma.
\]

## Cis-trans heterogeneity

Let \(w_c\) and \(w_t\) be normalized GIVW weight vectors for cis and trans subsets. Cross-block outcome covariance is

\[
\Omega_{ct}=D_{Y,c}R_{ct}D_{Y,t}.
\]

Then

\[
\mathrm{Cov}(\hat\theta_c,\hat\theta_t)
=w_c^T\Omega_{ct}w_t,
\]

\[
Z_{c-t}
=
\frac{\hat\theta_c-\hat\theta_t}
{\sqrt{SE_c^2+SE_t^2-2\mathrm{Cov}(\hat\theta_c,\hat\theta_t)}}.
\]

## Sibling co-perturbation omnibus test

For one sibling exposure with effect vector \(b\), standard-error matrix \(D\), and signed LD \(R\),

\[
\Sigma=DRD,
\]

\[
Q_{\mathrm{sib}}=b^T\Sigma^{-1}b,
\qquad
df=\mathrm{rank}(R).
\]

Sibling-test \(P\)-values are BH-adjusted across tested siblings.

## TF-locus diagnostics

For TF locus \(T\), leave-one-locus effect change is

\[
\Delta_T
=
\left|\hat\theta_{\mathrm{all}}-\hat\theta_{-T}\right|.
\]

Exposure-information weight:

\[
w_{T,\mathrm{exp}}
=
\frac{I_{\mathrm{all}}-I_{-T}}{I_{\mathrm{all}}}.
\]

MR-precision weight:

\[
w_{T,\mathrm{prec}}
=
\frac{J_{\mathrm{all}}-J_{-T}}{J_{\mathrm{all}}}.
\]

## Multivariable MR

For \(m\) retained SNPs and \(p\) exposures, let \(B\in\mathbb R^{m\times p}\) contain complete SNP-exposure effects and \(\Gamma_Y\in\mathbb R^m\) contain SNP-outcome effects.

\[
\Omega_Y=D_YRD_Y.
\]

Correlated-IV MV-IVW:

\[
\hat\theta
=
(B^T\Omega_Y^{-1}B)^{-1}
B^T\Omega_Y^{-1}\Gamma_Y.
\]

\[
\mathrm{Cov}(\hat\theta)
=
(B^T\Omega_Y^{-1}B)^{-1}.
\]

Residual statistic:

\[
Q_{MV}
=
(\Gamma_Y-B\hat\theta)^T
\Omega_Y^{-1}
(\Gamma_Y-B\hat\theta),
\qquad df=m-p.
\]

A standard MVMR fit requires every selected SNP-exposure pair in \(B\) and its standard-error matrix to be numerically observed.

## Scale-invariant MVMR condition diagnostic

With \(LL^T=\Omega_Y\),

\[
B_w=L^{-1}B.
\]

After column normalization,

\[
\widetilde B_{w,\cdot k}
=
\frac{B_{w,\cdot k}}{\|B_{w,\cdot k}\|_2},
\]

and the reported condition diagnostic is

\[
\kappa(\widetilde B_w).
\]

## Experimental correlated conditional-F diagnostic

For exposure \(i\), let \(\delta_i\) denote the fitted coefficients of exposure-association column \(i\) on the remaining columns. Define

\[
q_i=1,
\qquad
q_{-i}=-\delta_i.
\]

For SNP \(j\) and exposure \(k\),

\[
M_{jk}=SE_{jk}q_k.
\]

With exposure-estimation-error correlation matrix \(C\), the working residual covariance is

\[
V_i=(MCM^T)\circ R,
\]

where \(\circ\) is the Hadamard product. The residual vector is

\[
r_i=B_{\cdot i}-B_{\cdot,-i}\delta_i.
\]

The iterative diagnostic reports

\[
Q_i=r_i^TV_i^{-1}r_i,
\]

\[
F_{\mathrm{cond},i}^{\mathrm{experimental}}
=
\frac{Q_i}{m-p+1}.
\]

This statistic is labeled `EXPERIMENTAL_CORRELATED_IV_CONDITIONAL_F` in the implementation.

## References

1. Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC. Mendelian randomization with fine-mapped genetic data: Choosing from large numbers of correlated instrumental variables. *Genet Epidemiol.* 2017;41:714-725. https://doi.org/10.1002/gepi.22077
2. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *Int J Epidemiol.* 2019;48:713-727. https://doi.org/10.1093/ije/dyy262
