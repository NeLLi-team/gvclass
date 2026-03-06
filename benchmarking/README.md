# Benchmarking Workspace

This directory separates local benchmarking and reference-workspace assets from the production `resources/` runtime database.

Layout:

- `contamination/`
  - local contamination benchmark inputs, outputs, feature tables, and reports
- `completeness/`
  - local completeness/reference workdirs used to derive or evaluate completeness resources
- `checksums/`
  - optional private archive checksums
- `RUNBOOK.md`
  - exact commands and expectations for reproducing packaging and reruns

Policy:

- `resources/` stays production-focused and is the only directory intended for the end-user runtime download.
- `benchmarking/` is private/local for now and should not be published publicly until reviewed.
- If private transfer is needed, package `benchmarking/` separately and move it by `scp`, shared storage, or another internal channel.
