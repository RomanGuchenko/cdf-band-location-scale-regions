# Research Provenance and Attribution

## Project

**Title:** *From Distribution-Free CDF Bands to Polyhedral Confidence Regions for Location-Scale Families*  
**Artifact type:** Documented AI-generated research artifact  
**Human initiator, curator, and verifier:** Roman Guchenko  
**Substantive scientific development and initial manuscript generation:** ChatGPT (OpenAI, GPT-5.6 Sol)  
**Research developed:** August 2026  
**Initial public version:** 20 September 2026

## Purpose of this record

This document records how the research developed and how credit is assigned.

The project does not fit ordinary single-author terminology. The substantive mathematical and computational development was generated predominantly by ChatGPT during an extended dialogue with Roman Guchenko. The human participant initiated the investigation, selected among directions proposed during the dialogue, executed computational experiments, returned results to the model, later personally reviewed and verified the mathematics and R implementation, curated the resulting materials, and prepared the work for public dissemination.

The purpose of this provenance record is not to minimize either role. It is to make the division of labor explicit enough that a reader can evaluate the scientific artifact without having to infer who contributed what.

## 1. Origin of the investigation

The project began with a human question about whether confidence information for a normal cumulative distribution function, viewed relative to an empirical CDF, could be inverted to obtain confidence information for the normal location and scale parameters.

That initial question came from Roman Guchenko.

ChatGPT recognized the problem as one of inversion of distribution-free goodness-of-fit or CDF-band constraints. During the subsequent dialogue, the model developed the observation that suitable order-statistic constraints become affine inequalities in the location and scale parameters.

The model then generalized the construction beyond the initial normal-model setting to continuous two-parameter location-scale families and, later, beyond a single goodness-of-fit procedure to arbitrary simultaneous distribution-free order-statistic CDF bands.

The human participant did not specify those mathematical developments in advance. His role at this stage was to ask questions, choose among suggested research directions, request deeper investigation of particular possibilities, execute proposed computations, and return the numerical results to the model.

## 2. Development of the mathematical direction

The scientific direction changed during the conversation rather than following a predetermined plan.

Early in the investigation, ChatGPT treated some geometric properties of the inverted confidence region as potentially novel. Literature exploration then identified important precedents, including earlier work on goodness-of-fit inversion and KS-based location-scale confidence regions. The novelty claim was narrowed accordingly.

A later search identified recent work on inversion of distribution-free CDF bands through parametric quantile functions. This again changed the scope of the project.

The research then concentrated on the special geometric structure of two-parameter location-scale families:

- arbitrary order-statistic CDF bands produce affine constraints in `(mu, sigma)`;
- the induced confidence region is convex and polyhedral;
- the scale projection can be characterized through the gap between piecewise-affine envelopes;
- the location projection and area can be computed exactly;
- monotonicity of standard band endpoints permits an efficient envelope algorithm;
- the construction is location-scale equivariant;
- different CDF bands can be compared according to the parameter uncertainty they induce rather than according to CDF-band width alone.

This became the central technical direction of the paper.

The manuscript deliberately retains a cautious statement about novelty:

> A complete literature review is still required before any stronger novelty claim is made, especially for the general polyhedral and algorithmic statements.

## 3. Algorithms, code, and simulations

ChatGPT generated the principal derivations and drafts of the computational implementation during the dialogue.

The resulting R implementation includes:

- Monte Carlo calibration of Kolmogorov-Smirnov, Berk-Jones, and Duembgen-Wellner bands;
- analytic Dvoretzky-Kiefer-Wolfowitz bounds;
- transformation of order-statistic probability limits through a base quantile function;
- exact construction of the induced location-scale confidence region;
- envelope-based scale and location projections;
- exact joint area;
- handling of empty and unbounded regions;
- normalization for numerical stability and affine equivariance;
- Monte Carlo comparison of the four bands;
- reproduction drivers for the experiments reported in the paper.

The simulation study was designed in the dialogue and executed by the human participant using the generated implementation. Numerical outputs were returned to ChatGPT, which interpreted the results and incorporated them into the developing manuscript.

The reported experiments use 100,000 uniform calibration samples for the simulated critical values and 10,000 outer Monte Carlo replications. The paper compares the normal location-scale family across several sample sizes and also reports Laplace and Student-t2 experiments.

## 4. Error correction and iterative feedback

The preserved dialogue contains mistakes, corrections, and changes of direction.

One important example concerns the indexing of order-statistic CDF bounds. During development, an incorrect treatment of the empirical jump endpoints led to suspicious coverage behavior. The issue was traced to the distinction between the lower constraint based on `i/n` and the upper constraint based on `(i-1)/n` for a continuous CDF evaluated at an order statistic. The implementation and manuscript were corrected.

This episode is intentionally not removed from the provenance record. It documents an actual research feedback loop: a generated construction led to computation, the computation exposed a problem, and the subsequent reasoning corrected the construction.

More generally, the transcript preserves abandoned ideas and literature-driven narrowing of claims rather than presenting the final paper as though its argument and scope were known from the outset.

## 5. Manuscript generation

ChatGPT generated the initial manuscript text after the mathematical, computational, and simulation work had developed sufficiently.

The substantive AI contribution included, in broad terms:

- literature exploration;
- formulation and generalization of the framework;
- mathematical derivations and proofs;
- algorithm design;
- R code drafts;
- simulation design;
- interpretation of the numerical results;
- organization of the technical note;
- initial manuscript drafting.

The public manuscript has subsequently been curated for release and contains an explicit provenance statement before the technical paper.

The technical prose avoids ordinary joint-author language where it would imply conventional human authorship. The document instead uses formulations such as "this paper studies," "the study investigates," and "exact representations are derived."

## 6. Human contribution

Roman Guchenko's role is described as **human initiator, curator, and verifier**.

His contributions include:

- posing the initial research question;
- deciding which proposed directions were worth pursuing;
- requesting deeper development of selected ideas;
- executing computational experiments;
- returning numerical outputs into the dialogue;
- making decisions about the scope and continuation of the investigation;
- personally reviewing and verifying the mathematical arguments and R implementation before dissemination;
- preserving the transcript and computational materials;
- curating the final public research artifact;
- arranging publication and long-term preservation.

He does not claim conventional authorship of the scientific contributions that were generated by ChatGPT.

## 7. Contribution map

| Activity | Roman Guchenko | ChatGPT |
| --- | --- | --- |
| Initial research question | Primary | - |
| Selection of follow-up directions | Primary / interactive | Proposed alternatives and developments |
| Literature exploration | Reviewed as part of the project | Primary during the research dialogue |
| General mathematical formulation | Verification | Primary |
| Proofs and derivations | Verification | Primary |
| Algorithm design | Verification | Primary |
| R implementation drafts | Execution, testing, verification | Primary drafting |
| Simulation design | Selection, execution | Primary proposal and design |
| Simulation execution | Primary | - |
| Interpretation of results | Review and verification | Primary initial interpretation |
| Initial manuscript drafting | Review and curation | Primary |
| Final public curation | Primary | Assistance |
| Decision to disseminate | Primary | - |

The table is descriptive rather than a conventional authorship taxonomy.

## 8. Human verification

After the AI-generated scientific development, Roman Guchenko personally reviewed and verified the mathematical arguments and implementation before deciding to publish the work.

The repository additionally contains `location_scale_checks.R`, a separate verification suite. It is designed to test the implementation independently of the main production path.

The checks include:

1. envelope-evaluation regression tests, including one-line, final-segment, and nearly-parallel cases;
2. randomized comparison of the linear envelope scan with an independent quadratic pairwise formula;
3. direct envelope comparisons and dense numerical integration of the region area;
4. affine-equivariance checks under very small and very large changes of measurement units;
5. explicit empty and unbounded confidence-region cases;
6. equality between each inverted band event and its originating statistic event for KS, Berk-Jones, and Duembgen-Wellner constructions;
7. exchangeable-rank calibration, seeded reproducibility, and reuse of precomputed band calibrations;
8. an empirical log-log scaling diagnostic intended to catch an accidental return to quadratic geometry.

The automated suite is part of the reproducibility record. It is separate from the human curator's own responsibility for reviewing the mathematics and implementation.

## 9. Research transcript

The human-visible research dialogue is preserved because the provenance of the research is itself part of the public claim of this project.

Two transcript renderings are provided.

### Reader edition

`simulate-confidence-intervals-reader-edition.pdf`

This version renders mathematical notation and Markdown for readability while retaining the chronological sequence of the dialogue.

### Archival transcript

`simulate-confidence-intervals-transcript.pdf`

This version preserves the retrieved conversation record more directly.

The transcript is evidence for:

- the origin of the initial question;
- which directions were suggested by the model;
- which directions were selected by the human participant;
- mathematical derivations as they emerged;
- literature searches that narrowed earlier claims;
- mistakes and corrections;
- code development;
- simulation design and interpretation;
- the eventual decision to draft and publish the research.

The transcript contains the visible dialogue as retrieved from the source conversation. It does not purport to contain hidden model reasoning or internal chain-of-thought.

## 10. Known transcript limitation

The retrieved source conversation contains one known truncation.

One assistant message exceeded the source retrieval interface's 20,000-character limit. Both transcript renderings identify the point at which this occurred. The remainder of that message was not available through the retrieval path used to construct the transcript.

The missing text has not been guessed, regenerated, or silently reconstructed.

If a complete export of the original conversation becomes available, the preferred procedure is to restore the missing passage in a new archived version while retaining the earlier public version in version history.

Accordingly, references to the record should use terms such as **preserved record**, **available research record**, or **retrieved transcript**, rather than claiming that the present transcript is complete without qualification.

## 11. What this attribution does and does not claim

This project makes two distinctions explicit.

First, the research was **not autonomous in the sense of an AI system independently choosing a topic, operating a computer without a human interlocutor, and deciding on its own to publish the result**. The preserved dialogue shows continuing human interaction, direction selection, execution of experiments, and feedback.

Second, that human interaction is **not treated as evidence that the human participant originated the principal mathematical constructions, derivations, algorithms, code drafts, simulation design, interpretations, or initial manuscript**. Those contributions were generated predominantly by ChatGPT.

The chosen attribution is therefore:

> **Human initiator, curator, and verifier:** Roman Guchenko  
> **Substantive scientific development and initial manuscript generation:** ChatGPT (OpenAI, GPT-5.6 Sol)

This wording is intended to preserve both human accountability for curation and verification and accurate credit for the origin of the scientific development.

## 12. Public research record

The public research record is intended to include:

- the technical manuscript in PDF form;
- the LaTeX manuscript source;
- the reader-oriented transcript;
- the archival transcript;
- the R reference implementation;
- the independent R verification suite;
- the Jupyter notebook;
- the rendered notebook;
- the settings and seeds needed to reproduce the reported experiments;
- this provenance statement.

The full R reference implementation is also reproduced inside the PDF. This is intentional: the PDF is meant to remain a substantially self-contained scientific artifact rather than relying entirely on auxiliary repository files.

## 13. Versioning and preservation

The initial public artifact is intended to be released as:

**v1.0 - 20 September 2026**

Subsequent corrections, recovered transcript material, expanded simulations, or changes to the scientific manuscript should be issued as new versions rather than silently replacing the historical record.

The Git repository is intended to preserve development and source history. Tagged releases should identify frozen public versions of the artifact.
