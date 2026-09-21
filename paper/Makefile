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
#   make tables   regenerate paper/tables/ from config.yml and the data
#
# refs.bib is generated from the "postdoc" collection of the Zotero library.
# Add or correct a reference in Zotero, then `make refs` -- never edit the .bib.
# Reading Zotero's database takes no lock, so Zotero can stay open.
#
# If Better BibTeX ever gets installed, exporting that collection as Better
# BibLaTeX with "keep updated" ticked makes the file rewrite itself on every
# edit and `make refs` redundant. The citation keys match either way.
#
# The workflow figure is a hand-drawn SVG exported to PNG with Inkscape. Data
# figures, when the list is settled, go here as one rule per figure, depending
# on ../data/predictions.parquet and run under ~/.venvs/streamlit/bin/python
# (which has plotly and kaleido; the repo's .venv has neither).

ZOTERO := ../scripts/zotero_to_bib.py
PYTHON := ../.venv/bin/python
INKSCAPE := /Applications/Inkscape.app/Contents/MacOS/inkscape

paper.pdf: paper.tex refs.bib tables/parasites.tex figures/workflow.png
	tectonic --synctex --keep-logs paper.tex

figures/workflow.png: figures/workflow.svg
	$(INKSCAPE) $< --export-type=png --export-dpi=300 --export-background=white --export-filename=$@

# The filter funnel: proteins left after each stage of the pipeline, per species.
# Rereads the TISSUES downloads and the EggNOG/COG files, so it takes a minute.
figures/funnel.svg figures/funnel.pdf tables/funnel.csv: ../config.yml ../data/eligible_proteins.parquet ../data/predictions.parquet ../scripts/build_filter_funnel.py
	cd .. && $(PYTHON) scripts/build_filter_funnel.py

figures/funnel.png: figures/funnel.svg
	$(INKSCAPE) $< --export-type=png --export-dpi=300 --export-background=white --export-filename=$@

# The multi-host comparison: links shared between and specific to the hosts of each
# parasite with several, and why a host lacks a link another host carries.
figures/multi_host.svg figures/multi_host.pdf tables/multi_host.csv: ../config.yml ../data/predictions.parquet ../data/host_orthologs.parquet ../scripts/build_multi_host_figure.py
	cd .. && $(PYTHON) scripts/build_multi_host_figure.py

figures/multi_host.png: figures/multi_host.svg
	$(INKSCAPE) $< --export-type=png --export-dpi=300 --export-background=white --export-filename=$@

# The host families the most parasites converge on, over the parasites of human.
figures/shared_families.svg figures/shared_families.pdf tables/shared_families.csv: ../config.yml ../data/predictions.parquet ../scripts/build_shared_families.py
	cd .. && $(PYTHON) scripts/build_shared_families.py

figures/shared_families.png: figures/shared_families.svg
	$(INKSCAPE) $< --export-type=png --export-dpi=300 --export-background=white --export-filename=$@

figures: figures/workflow.png figures/funnel.pdf figures/funnel.png figures/multi_host.pdf figures/multi_host.png figures/shared_families.pdf figures/shared_families.png

# The parasite table, in three renderings: .tex for the paper, .csv and .html
# for pasting into Google Sheets / Docs.
tables/parasites.tex: ../config.yml ../data/predictions.parquet ../scripts/build_parasite_table.py
	cd .. && $(PYTHON) scripts/build_parasite_table.py

tables: tables/parasites.tex

# Not a dependency of paper.pdf: the Zotero library changes under the Makefile's
# feet, and a build that silently rewrote the bibliography would be worse than
# one that used yesterday's.
refs:
	../.venv/bin/python $(ZOTERO) --collection postdoc --out refs.bib

clean:
	rm -f paper.pdf paper.aux paper.bbl paper.bcf paper.blg paper.log \
	      paper.out paper.run.xml paper.synctex.gz

.PHONY: clean refs tables figures
