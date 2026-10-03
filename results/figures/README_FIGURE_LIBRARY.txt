MIDUS TWO-PHASE SIMULATION — FIGURE LIBRARY

Key changes:
  - relative precision now uses one common -25% to 450% scale across all scenarios
  - relative precision includes a dotted zero reference line
  - enlarged white mask behind AIPW so connecting lines do not show through
  - common power scale changed to 20%-100% with the 80% reference retained
  - runtime method order explicitly locked to Naive, CCA, FCS-MI, JM-MI, IPW, AIPW
  - Phase-2 sampling fraction displayed with phi_2 to avoid conflict with selection-model eta notation
  - AIPW uses native shape 10 (open circle-plus) in both panels and legend
  - null rejection-rate axis explicitly labels 4%, 5%, and 6%
  - EmpSE and ModSE abbreviations standardized in SE figures
  - runtime methods displayed in the same top-to-bottom order as the manuscript figures
  - Dark2 method colors
  - five white-filled open symbols plus native plotting symbol 10 for AIPW
  - thicker method lines and darker text for manuscript-ready output
  - one-row method legend
  - equally spaced Phase-1 sample-size design levels
  - two-line Phase-2 facet strips using phi_2
  - PDF-safe plotmath theta and phi labels
  - cairo_pdf output
  - fixed y-axis scales by performance measure (bias: -0.01 to 0.025)
  - no shaded acceptable-performance bands
  - reference boundary lines for SE ratio, coverage, and Type I error
  - redesigned horizontal runtime plot with explicit log-scale tick labels

01_full_factorial/
  Complete R_A^2 x R_Y^2 figure library.

02_manuscript_candidates/
  Weak, MIDUS-like, and strong auxiliary-information candidate figures.

03_final_manuscript/
  Parsimonious manuscript-facing selection of core, stress, and supporting figures.

04_diagnostics/
  Redesigned runtime figure; failure rates are stored as data only because all were zero.

results/figure_data/
  Data underlying figures plus a figure manifest.
