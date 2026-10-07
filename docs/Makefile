# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# To run all notebooks that generate files like csv that are needed to populate rst files in docs
# NOTEBOOKS := $(filter-out jupyter_notebooks/MagmaSat_Crashing.ipynb,$(wildcard jupyter_notebooks/*.ipynb))

# .PHONY: notebooks
# notebooks:
# 	@failed=""; \
# 	for nb in $(NOTEBOOKS); do \
# 		echo "Running $$nb"; \
# 		jupyter nbconvert --execute --to notebook --inplace "$$nb" || failed="$$failed $$nb"; \
# 	done; \
# 	git --no-pager diff --stat -- tables; \
# 	if [ -n "$$failed" ]; then echo "FAILED:$$failed"; exit 1; fi
# 	@$(MAKE) html

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
