# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Custom targets for local development
.PHONY: serve livehtml clean-serve

serve: html
	@echo "Starting local server at http://localhost:8000"
	@echo "Press Ctrl+C to stop"
	cd $(BUILDDIR)/html && python -m http.server 8000

livehtml:
	@echo "Starting live reload server (requires sphinx-autobuild)"
	@echo "Install with: pip install sphinx-autobuild"
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" --open-browser --port 8000

clean-serve: clean serve

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
