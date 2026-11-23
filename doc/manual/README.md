# FoamWaveTank Manual

`FoamWaveTank_Manual.pdf` is the compiled manual. The LaTeX source is
`FoamWaveTank_Manual.tex` (front matter), `chapters/*.tex` and `references.bib`.

The manual covers only what FoamWaveTank adds to waves2Foam; for everything
else it refers to the waves2Foam Manual (N. G. Jacobsen, 2017).

## Building the PDF

Needs a LaTeX distribution (on Ubuntu: `texlive-latex-extra`, `texlive-fonts-recommended`,
`latexmk`), or upload the folder to Overleaf.

```bash
cd doc/manual
latexmk -pdf FoamWaveTank_Manual.tex
latexmk -c          # remove auxiliary files
```

Commit the updated `FoamWaveTank_Manual.pdf` together with the source changes.
