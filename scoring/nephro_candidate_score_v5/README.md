# Nephro Candidate Score v5

## Allele Frequency Substitution & Model Approximation

### Background
In `scoring/nephro_candidate_score/` (v4 model), the variant scoring formula separates exome and genome allele frequencies:
```text
((gnomade_variant - 0.00658) / 0.05959) * (-309.33539) + ((gnomadg_variant - 0.02425) / 0.11003) * (-2.54581)
```
where:
- `gnomade_variant` ($E$) represents dbNSFP exome allele frequency
- `gnomadg_variant` ($G$) represents dbNSFP genome allele frequency

### v5 Approximation ($E = G = J$)
In `scoring/nephro_candidate_score_v5/`, gnomAD v4.1 data provides a joint frequency:
- `dbNSFP_gnomAD4.1_joint_AF` mapped to `gnomad_joint_variant` ($J$)

The v5 formula substitutes `gnomad_joint_variant` for both terms:
```text
((gnomad_joint_variant - 0.00658) / 0.05959) * (-309.33539) + ((gnomad_joint_variant - 0.02425) / 0.11003) * (-2.54581)
```
This produces one effective linear contribution:
```text
-5214.199508 * J + 34.718271
```

### Rationale: Why Terms Must Not Be Mechanically Deleted
Mechanically deleting either AF term changes both the effective slope and intercept of the model:
- Deleting the exome term leaves only the genome coefficient (slope -23.137, intercept 0.561).
- Deleting the genome term leaves only the exome coefficient (slope -5191.062, intercept 34.157).

Preserving the intended response under equal AF inputs ($E = G = J$) requires maintaining the joint approximation until a replacement model is re-fitted and statistically validated against clinical benchmarks.

### Verified Reference Outputs (Regression Contract)
1. **Equal-AF Case ($E = G = J = 0.0005$):**
   - Inputs: CADD = 20.0, EFFECT = missense_variant, IMPACT = MODERATE, ngs = 0.0, pattern = unknown.
   - Legacy score: `0.056695344`
   - v5 score: `0.056695344` (identical to legacy)

2. **Unequal-AF Case ($E = 0.0, G = 0.001, J = 0.0005$):**
   - Inputs: Same as above.
   - Legacy score: `0.443311077`
   - v5 score: `0.056695344`
