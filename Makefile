# Agents vs humans (TREM2 paper) — Makefile (mirror of mise.toml for `make` users)
#
# Usage:
#   make setup            # install deps
#   make analysis-all     # run every canonical analysis
#   make figures-blog     # serve the 7 blog HTML figures locally
#   make help

PY ?= uv run python

.PHONY: help setup lint format \
        analysis-human-vs-agent analysis-leaderboard analysis-diversity analysis-methods \
        analysis-all figures-blog figures-paper clean

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk -F':.*?## ' '{printf "  %-28s %s\n", $$1, $$2}'

setup: ## Install Python deps via uv
	uv sync

lint: ## Lint with ruff
	uv run ruff check .

format: ## Auto-format with ruff
	uv run ruff format .

analysis-human-vs-agent: ## Headline cohort stats
	$(PY) analyses/human_vs_agent/headline_stats.py

analysis-leaderboard: ## Per-design top-N leaderboards
	$(PY) analyses/leaderboard/build_leaderboards.py

analysis-diversity: ## Per-team / cohort sequence diversity
	$(PY) analyses/sequence_diversity/diversity.py

analysis-methods: ## Design method × outcome cross-tabs
	$(PY) analyses/methods/method_outcome_xtab.py

analysis-all: analysis-human-vs-agent analysis-leaderboard analysis-diversity analysis-methods ## Run every canonical analysis

figures-blog: ## Serve the 7 hand-authored blog HTML figures on http://localhost:8080
	python -m http.server -d figures/blog 8080

figures-paper: ## Render matplotlib paper figures (placeholder — populate scripts/plotting/)
	@echo "No matplotlib figure pipeline wired up yet. Add scripts to scripts/plotting/ and edit this rule."
	@echo "The 7 blog HTML figures already live in figures/blog/ — see 'make figures-blog'."

clean: ## Remove derived outputs (does NOT touch data/)
	find analyses -name "report.md" -delete
	find analyses -name "*.json" -delete
	find analyses -name "top*.csv" -o -name "*_winners.csv" -delete
