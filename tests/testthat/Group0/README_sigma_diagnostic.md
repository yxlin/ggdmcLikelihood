# Sigma Likelihood Profile Diagnostic

## Purpose

This diagnostic script (`05_sigma_diagnostic.r`) tests whether the CDM likelihood function with MVN attribute distribution correctly identifies the true `sigma` parameter value.

The user hypothesis being tested is: **When all other parameters are fixed at their true values and only sigma varies, the likelihood might be larger at non-true sigma values**, which would indicate a bug in the likelihood implementation.

## What the Script Does

1. **Generates synthetic data** with known true parameter values:
   - True sigma = 0.2
   - True means = [0.5, 0.2]
   - True guess and slip parameters (varying by item)

2. **Fixes all parameters at true values** except sigma

3. **Evaluates likelihood across a grid** of sigma values (0.01 to 0.95 in steps of 0.02)

4. **Tests multiple sample sizes**: N = 100, 500, 2000, 5000
   - Larger N should converge better to true parameter if likelihood is correct

5. **Creates diagnostic plots** showing likelihood profile vs sigma

6. **Reports test results**:
   - PASS: Maximum likelihood occurs at or very near true sigma
   - FAIL: Maximum likelihood occurs at a different sigma value (suggests bug)
   - CHECK: Marginal case requiring investigation

## How to Run

```bash
cd /media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat
Rscript Group0/05_sigma_diagnostic.r
```

## Output Files

- **Console output**: Detailed results for each sample size
- **PDF plot**: `Group0/sigma_likelihood_profile.pdf` - Shows likelihood curves for all sample sizes

## Interpreting Results

### Expected Behavior (No Bug)

If the likelihood function is correctly specified:
- The maximum log-likelihood should occur **at or very near** the true sigma (0.2)
- As N increases, the maximum should converge closer to true sigma
- The likelihood curve should be smooth and peaked at the true value

### Warning Signs (Potential Bug)

If you see:
- Maximum likelihood at a **different sigma value** (especially if consistent across sample sizes)
- **Flat or multi-modal** likelihood profiles
- **Non-monotonic convergence** as N increases (larger N getting worse fit)
- **Consistent bias** away from true value

These would suggest the MVN profile probability calculation has an issue.

## Test Criteria

The script classifies each test as:

- **PASS**: `|LogLik_Diff| < 0.01` - Maximum is at/near true sigma
- **FAIL**: `LogLik_Diff > 0.01` - Maximum is at a different sigma (suggests bug)
- **CHECK**: `LogLik_Diff < -0.01` - Marginal case (may need investigation)

Where `LogLik_Diff = LogLik_at_ML_sigma - LogLik_at_true_sigma`

## Additional Diagnostics

The script also reports:
1. **Curvature at true sigma**: Checks if true sigma is a local maximum
2. **Likelihood values in neighborhood**: Shows nearby sigma values and their likelihoods
3. **Summary table**: Compact view of all results

## What to Do if Tests Fail

If multiple tests show FAIL status:

1. Check if the pattern is consistent (same direction of error)
2. Examine the plot - is there a clear maximum elsewhere?
3. Check if larger N makes it worse (suggests systematic bias)
4. Investigate the `calculate_mvn_profile_probabilities()` function in `cdm.h`
5. Verify the pmvnorm() calls are using correct bounds
6. Check if the compound symmetry correlation matrix is correctly constructed

## Notes

- The script uses fixed seed (123 + i) for reproducibility
- All parameters except sigma are fixed at true values during profiling
- Uses the same DINO rule and 5-item Q-matrix as the original test
- Tests use the caching optimization recently added to `cdm.h`
