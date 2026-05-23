# OSCAR — Product Context

## Product Purpose
Internal lab tooling for the Romagnani lab at Charité. Helps researchers prepare single-cell sequencing experiments: generating metadata files and antibody reference CSVs before running the OSCAR Nextflow pipeline.

## Users
Biology PhD students and postdocs at Charité. Primarily wet lab researchers who understand their experiments but are not bioinformaticians. They know their assay type, cell populations, and antibody panels — but may not know what "dual index" or "TotalSeq-B" means in a pipeline context. Need enough explanation to fill in fields correctly without needing to ask a bioinformatician.

## Register
product

## Tone
Clear, calm, professional. No jargon without explanation. No patronising over-explanation. A student should be able to fill in the form correctly on their first attempt if they know their own experiment.

## Anti-references
- SaaS-dashboard chrome (sidebar navs, metrics panels, activity feeds)
- Clinical/EHR aesthetic (white + teal, form-heavy corporate)
- Developer-tool dark mode

## Design principles
- The tool should disappear into the task. A student should finish in under 5 minutes.
- Every field label and hint must answer "what do I put here?" not "what is this concept?"
- Errors should tell you how to fix them, not just what went wrong.
- The download button must always be reachable without scrolling.
