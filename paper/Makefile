# OrthoHPI 2.0 manuscript.
#
#   make          build paper.pdf
#   make clean    remove build products
#
# Needs tectonic (brew install tectonic). It fetches the LaTeX packages the
# document uses on first build -- that run is slow, later ones are not -- and
# runs the multi-pass and biber loop itself, so there is nothing else to call.
#
#   make refs     re-export refs.bib from Zotero
#
# refs.bib is generated from the "postdoc" collection of the Zotero library.
# Add or correct a reference in Zotero, then `make refs` -- never edit the .bib.
# Reading Zotero's database takes no lock, so Zotero can stay open.
#
# If Better BibTeX ever gets installed, exporting that collection as Better
# BibLaTeX with "keep updated" ticked makes the file rewrite itself on every
# edit and `make refs` redundant. The citation keys match either way.
#
# There are no figure rules yet. When the figure list is settled they go here
# as one rule per figure, depending on ../data/predictions.parquet and run
# under ~/.venvs/streamlit/bin/python (which has plotly and kaleido; the repo's
# .venv has neither).

ZOTERO := ../scripts/zotero_to_bib.py

paper.pdf: paper.tex refs.bib
	tectonic --synctex --keep-logs paper.tex

# Not a dependency of paper.pdf: the Zotero library changes under the Makefile's
# feet, and a build that silently rewrote the bibliography would be worse than
# one that used yesterday's.
refs:
	../.venv/bin/python $(ZOTERO) --collection postdoc --out refs.bib

clean:
	rm -f paper.pdf paper.aux paper.bbl paper.bcf paper.blg paper.log \
	      paper.out paper.run.xml paper.synctex.gz

.PHONY: clean refs
