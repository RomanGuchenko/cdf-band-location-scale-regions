# From Distribution-Free CDF Bands to Polyhedral Confidence Regions for Location-Scale Families

**A documented AI-generated research artifact**

**Human initiator, curator, and verifier:** Roman Guchenko  
**Substantive scientific development and initial manuscript generation:** ChatGPT (OpenAI, GPT-5.6 Sol)

**Public version:** 20 September 2026

## Overview

This repository preserves the public research record for the note *From Distribution-Free CDF Bands to Polyhedral Confidence Regions for Location-Scale Families*.

The paper studies confidence regions for the location and scale parameters of continuous location-scale families obtained by inverting simultaneous, distribution-free CDF bands. The resulting order-statistic constraints become affine inequalities in the location and scale parameters, giving convex polyhedral confidence regions. The note develops exact formulas and algorithms for the scale and location projections, piecewise-linear boundaries, region area, and an efficient monotone-slope envelope construction. It also compares Kolmogorov-Smirnov, Dvoretzky-Kiefer-Wolfowitz, Berk-Jones/Owen, and Duembgen-Wellner bands by the parametric uncertainty they induce.

The unusual feature of this project is its provenance. The substantive mathematical and computational development was generated predominantly by ChatGPT during an extended research dialogue. Roman Guchenko supplied the initial research question, selected directions for further investigation, executed computational experiments, returned results to the model, and later personally reviewed and verified the mathematical arguments and R implementation. He curates and publishes the resulting research record without claiming conventional authorship of the AI-generated scientific contributions.

For the full account, see [`PROVENANCE.md`](PROVENANCE.md).

## Repository contents

| Path | Description |
| --- | --- |
| [`ls_note2_revised.pdf`](paper/ls_note2_revised.pdf) | Public PDF of the research note |
| [`ls_note2_revised.tex`](paper/ls_note2_revised.tex) | LaTeX source |
| [`simulate-confidence-intervals-reader-edition.pdf`](transcript/simulate-confidence-intervals-reader-edition.pdf) | Reader-oriented rendering of the chronological human-AI dialogue |
| [`simulate-confidence-intervals-transcript.pdf`](transcript/simulate-confidence-intervals-transcript.pdf) | Archival transcript as retrieved from the source conversation |
| [`location_scale_reference.R`](code/location_scale_reference.R) | Pedagogical reference implementation and reproduction drivers |
| [`location_scale_checks.R`](code/location_scale_checks.R) | Independent verification and regression checks |
| [`location_scale_research.ipynb`](notebook/location_scale_research.ipynb) | Computational notebook |
| [`location_scale_research.html`](notebook/location_scale_research.html) | Rendered notebook |
| [`PROVENANCE.md`](PROVENANCE.md) | Detailed research provenance and attribution record |

The R reference implementation is also reproduced in the PDF so that the main research artifact remains computationally inspectable even when read independently of this repository.

## Main scientific contents

For a continuous location-scale family,

```text
F_{mu,sigma}(x) = G((x - mu)/sigma), sigma > 0,
```

simultaneous order-statistic bounds

```text
ell_i <= F(X_(i)) <= u_i
```

imply affine constraints

```text
X_(i) - sigma G^{-1}(u_i) <= mu <= X_(i) - sigma G^{-1}(ell_i).
```

The paper uses this representation to study:

- convex and polyhedral confidence regions in `(mu, sigma)`;
- exact scale and location projections;
- piecewise-linear boundaries and exact joint area;
- a linear-time envelope construction after sorting when band endpoints are monotone;
- location-scale equivariance;
- Monte Carlo comparisons of several distribution-free CDF bands;
- the broader problem of designing distribution-free bands for parameter-specific objectives.

The manuscript deliberately limits its novelty claims. In particular, it notes that goodness-of-fit inversion and KS-based location-scale confidence regions have important precedents, and that a complete literature review is still required before any stronger novelty claim is made for the general polyhedral or algorithmic statements.

## Verification

The human curator personally reviewed and verified the mathematics and implementation before public dissemination.

The repository also contains a separate automated verification suite. Among other checks, it:

- compares the envelope implementation with an independent quadratic pairwise formula on randomized inputs;
- compares envelope values with direct maxima and minima;
- checks exact area against dense numerical integration;
- tests affine equivariance under extreme changes of units;
- exercises explicit empty and unbounded confidence-region cases;
- checks that inversion reproduces the originating band-statistic event;
- checks rank calibration, seeded reproducibility, and reuse of precomputed bands;
- reports an empirical log-log timing slope as a guard against accidental quadratic behavior.

To run the verification script, use the `code/` directory as the working directory and run:

```bash
Rscript location_scale_checks.R
```

## Reproduction

The reference implementation contains the complete computational pipeline used in the note. The function `run_reproduction_study()` reproduces the reported simulation design: the required distribution-free bands are calibrated once for each sample size and then reused across the normal, Laplace, and Student-t2 location-scale families.

The default reproduction settings in the source use 100,000 calibration samples for simulated critical values and 10,000 outer Monte Carlo replications, so a complete run may be computationally substantial.

## Transcript and provenance record

The research-development transcript is included because provenance is part of the scientific claim of this artifact, not merely supplementary commentary.

Two versions are preserved:

- the **reader edition**, which renders mathematical notation and Markdown while retaining the chronological dialogue;
- the **archival transcript**, which preserves the retrieved source record more directly.

The retrieved source record contains one known truncation: one assistant message exceeded the source retrieval interface's 20,000-character limit. The missing continuation has not been silently reconstructed. If the original conversation export becomes available, the archival record can be updated with the recovered passage while preserving version history.

## Attribution

This repository uses role labels rather than conventional single-author attribution:

> **Human initiator, curator, and verifier:** Roman Guchenko  
> **Substantive scientific development and initial manuscript generation:** ChatGPT (OpenAI, GPT-5.6 Sol)

This wording is intended to describe the research process accurately. The work was not produced autonomously without human interaction, but the human role should also not be interpreted as a claim that the principal mathematical constructions, derivations, algorithms, code drafts, simulation design, or initial manuscript were originated by the human participant.

See [`PROVENANCE.md`](PROVENANCE.md) for the detailed record.

## Versioning

The initial public research artifact is intended to be tagged as:

**v1.0 - 20 September 2026**

Later corrections or additions should be released as new tagged versions rather than silently replacing the historical record.

## Licensing

This repository uses **mixed licensing** because it contains both software and
research/documentary material.

- **Software and code** are released under the **MIT License**, to the extent
  that licensable rights subsist and are held by Roman Guchenko. See
  [`LICENSE-CODE`](LICENSE-CODE).
- **The manuscript, repository documentation, provenance material, and other
  non-software research materials** are made available under
  **Creative Commons Attribution 4.0 International (CC BY 4.0)**, to the
  extent that licensable rights subsist and are held by Roman Guchenko. See
  [`LICENSE-DOCUMENTS`](LICENSE-DOCUMENTS).

The licensing notices do **not** assert that copyright necessarily subsists
in AI-generated material, or that Roman Guchenko owns exclusive copyright in
such material. They also do not override third-party rights in quoted,
referenced, or otherwise separately protected material.

Suggested attribution for the research artifact:

> *From Distribution-Free CDF Bands to Polyhedral Confidence Regions for Location-Scale Families*  
> Roman Guchenko — human initiator, curator, and verifier  
> Scientific development and initial manuscript generation: ChatGPT (OpenAI, GPT-5.6 Sol)  
> 2026
