## System details:

#### Positron and OS details:

```
Positron Version: 2026.03.0 build 212
Code - OSS Version: 1.108.0
Commit: f3aae65e0a1a11d39226cd884520f49301daef82
Date: 2026-02-26T05:01:47.944Z
Electron: 39.2.7
Chromium: 142.0.7444.235
Node.js: 22.21.1
V8: 14.2.231.21-electron.0
OS: Linux x64 6.8.0-106-generic
```

#### Session details:

```
R 4.5.1 

Posit.air-vscode Version 0.20.0
```

## Describe the issue:

The "Format Selection" command (`Ctrl+Shift+A`) formats the entire code cell in a Quarto document (`.Qmd`) rather than only the highlighted/selected code. This occurs in R code chunks — regardless of how much text is selected, the formatter applies to the full cell.

## Steps to reproduce the issue:

<!-- Reproduced in a Quarto (.Qmd) document with R code chunks -->

1. Open a `.Qmd` file containing R code chunks.
2. Highlight a subset of code within an R code chunk (e.g., a single line or a few lines, but not the entire chunk).
3. Run `Format Selection` via `Ctrl+Shift+A` (or via the Command Palette).
4. Observe that the entire code chunk is reformatted, not just the selected portion.

## Expected or desired behavior:

Only the highlighted/selected code should be reformatted. Code outside the selection within the same chunk should remain untouched. This is consistent with how "Format Selection" works in VS Code for other languages.

## Were there any error messages in the UI, Output panel, or Developer Tools console?

No error messages. The command executes without errors — it just applies to a broader scope than expected.

---

*This issue was drafted with the assistance of Claude Opus 4.6.*
