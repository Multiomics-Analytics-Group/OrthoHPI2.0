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
# The workflow figure is a hand-drawn SVG exported to PNG with Inkscape. Every
# other figure is one rule per build script, which writes the 300-dpi .png the
# paper includes and, in the same run, an .svg of it.

ZOTERO := ../scripts/zotero_to_bib.py
# the recipes below run from the repo root (`cd ..`), so this is relative to it
PYTHON := .venv/bin/python
INKSCAPE := /Applications/Inkscape.app/Contents/MacOS/inkscape

# The figures the manuscript includes, all as PNGs.
FIGURES := figures/workflow.png figures/funnel.png figures/networks.png \
	figures/parasite_similarity.png figures/shared_families.png \
	figures/multi_host.png figures/case_study.png figures/tissues_coverage.png \
	figures/confidence.png

paper.pdf: paper.tex refs.bib tables/parasites.tex $(FIGURES)
	tectonic --synctex --keep-logs paper.tex

figures/workflow.png: figures/workflow.svg
	$(INKSCAPE) $< --export-type=png --export-dpi=300 --export-background=white --export-filename=$@

# The filter funnel: proteins left after each stage of the pipeline on the host and
# parasite side of each host-parasite pair. Rereads the TISSUES downloads and the
# EggNOG members file, so it takes a minute.
figures/funnel.svg figures/funnel.png tables/funnel.csv tables/funnel_supplementary.csv: ../config.yml ../data/eligible_proteins.parquet ../data/predictions.parquet ../scripts/build_filter_funnel.py
	cd .. && $(PYTHON) scripts/build_filter_funnel.py

# The multi-host comparison: links shared between and specific to the hosts of each
# parasite with several, and why a host lacks a link another host carries.
figures/multi_host.svg figures/multi_host.png tables/multi_host.csv: ../config.yml ../data/predictions.parquet ../data/host_orthologs.parquet ../scripts/build_multi_host_figure.py
	cd .. && $(PYTHON) scripts/build_multi_host_figure.py

# The host families the most parasites converge on, over the parasites of human.
figures/shared_families.svg figures/shared_families.png tables/shared_families.csv: ../config.yml ../data/predictions.parquet ../scripts/build_shared_families.py
	cd .. && $(PYTHON) scripts/build_shared_families.py

# How alike the interactomes of the parasites of human are, pair by pair.
figures/parasite_similarity.svg figures/parasite_similarity.png tables/parasite_similarity.csv tables/parasite_similarity_medians.csv: ../config.yml ../data/predictions.parquet ../scripts/build_parasite_similarity.py
	cd .. && $(PYTHON) scripts/build_parasite_similarity.py

# The case study: one host-parasite pair from parasite family to host cell type.
figures/case_study.svg figures/case_study.png tables/case_study_links.csv tables/case_study_cell_types.csv: ../config.yml ../data/predictions.parquet ../data/tissues_cell_types.parquet ../scripts/build_case_study.py
	cd .. && $(PYTHON) scripts/build_case_study.py

# The predicted networks at a glance: interactions per parasite by host and by
# host-protein location.
figures/networks.svg figures/networks.png tables/networks.csv: ../config.yml ../data/predictions.parquet ../data/deeploc_localisations.parquet ../scripts/build_networks_figure.py
	cd .. && $(PYTHON) scripts/build_networks_figure.py

# Supplementary: how many STRING proteins of each host TISSUES annotates to each tissue
# at score 1. Rereads the TISSUES downloads.
figures/tissues_coverage.svg figures/tissues_coverage.png tables/tissues_coverage.csv: ../config.yml ../scripts/build_tissues_coverage.py
	cd .. && $(PYTHON) scripts/build_tissues_coverage.py

# Supplementary: the confidence scores of the predictions of every host-parasite pair.
figures/confidence.svg figures/confidence.png tables/confidence.csv: ../config.yml ../data/predictions.parquet ../scripts/build_confidence_figure.py
	cd .. && $(PYTHON) scripts/build_confidence_figure.py

figures: figures/workflow.png figures/funnel.png figures/multi_host.png figures/shared_families.png figures/parasite_similarity.png figures/case_study.png figures/networks.png figures/tissues_coverage.png figures/confidence.png

# The parasite table, in three renderings: .tex for the paper, .csv and .html
# for pasting into Google Sheets / Docs.
tables/parasites.tex: ../config.yml ../scripts/build_parasite_table.py
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
