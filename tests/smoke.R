source("CATER_MR.R")

# Column standardization and p-value derivation.
ss <- data.frame(
  SNP = c("rs1", "rs2"), CHR = c(1, 1), BP = c(100, 200),
  A1 = c("A", "C"), A2 = c("G", "T"), EAF = c(.2, .3),
  BETA = c(.2, -.1), SE = c(.05, .04), N = c(500, 500)
)
s <- .cater_standardize_sumstats(ss)
stopifnot(nrow(s) == 2L, all(is.finite(s$p)), identical(s$snp, c("rs1", "rs2")))

# One-hop region map.
grn <- data.frame(TF = c("TF1", "TF2"), Target = c("X", "X"))
ann <- data.frame(symbol = c("X", "TF1", "TF2"), chr = c("1", "1", "2"), tss = c(1000, 5000, 10000))
regions <- .cater_make_regions("X", c("TF1", "TF2"), ann, cis_window = 100, tf_window = 100)
q <- data.frame(snp = c("c", "t1", "t2"), chr = c("1", "1", "2"), pos = c(1000, 5000, 10000),
                a1 = "A", a2 = "G", beta = .1, se = .02, p = 1e-9, eaf = .2, n = 500)
cm <- .cater_candidate_map(q, regions, "X")
stopifnot(identical(cm$source, c("cis", "trans", "trans")))

# Allele harmonization: swapped outcome alleles flip beta.
exp <- q[1:2, ]
exp$source <- c("cis", "trans")
exp$parent_tf <- c("", "TF1")
out <- data.frame(SNP = c("c", "t1"), A1 = c("A", "G"), A2 = c("G", "A"),
                  BETA = c(.05, .03), SE = c(.01, .01))
out <- .cater_standardize_sumstats(out, label = "outcome", require_position = FALSE)
h <- .cater_harmonize(exp, out)
stopifnot(nrow(h) == 2L, abs(h$by[1] - .05) < 1e-12, abs(h$by[2] + .03) < 1e-12)

# Identity-LD GIVW matches ordinary fixed-effect IVW formula.
d <- data.frame(snp = c("a", "b"), source = "cis", parent_tf = "",
                bx = c(.2, .1), bx_se = c(.02, .02), by = c(.1, .06), by_se = c(.02, .03))
fit <- .cater_givw(d, diag(2))
w <- 1 / d$by_se^2
manual <- sum(w * d$bx * d$by) / sum(w * d$bx^2)
stopifnot(abs(fit$beta - manual) < 1e-12)

# Manc-COJO .ldr.cojo can contain chromosome blocks; cross-chromosome LD is zero.
f <- tempfile(fileext = ".ldr.cojo")
writeLines(c(
  "# Chromosome 1",
  "SNP\trs1\trs2",
  "rs1\t1\t0.25",
  "rs2\t0.25\t1",
  "# Chromosome 2",
  "SNP\trs3",
  "rs3\t1"
), f)
ld <- .cater_read_manc_ldr(f, c("rs1", "rs2", "rs3"))
stopifnot(abs(ld["rs1", "rs2"] - 0.25) < 1e-12,
          ld["rs1", "rs3"] == 0,
          ld["rs3", "rs3"] == 1)
unlink(f)

cat("CATER-MR smoke tests passed\n")
