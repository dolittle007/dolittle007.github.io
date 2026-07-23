# dolittle007.github.io

Personal website of **Ting-You Wang** — Bioinformatics Research Associate at Northwestern University Feinberg School of Medicine, working on cancer genomics, immunogenomics, and computational biology.

🔗 Live at [databeauty.com](https://databeauty.com)

## Stack

Built with [Jekyll](https://jekyllrb.com/) on the [al-folio](https://github.com/alshedivat/al-folio) v1 starter, which loads its runtime (layouts, search, comments, math, charts, CV rendering, etc.) from a set of pluggable `al_folio_*` gems rather than owning that code directly. See [AGENTS.md](AGENTS.md) and [docs/BOUNDARIES.md](docs/BOUNDARIES.md) for how starter vs. plugin ownership is split.

## Local development

```bash
bundle install
bundle exec jekyll serve
```

Visit `http://localhost:4000/`.

To build a production copy:

```bash
bundle exec jekyll build
```

## Content

- `_posts/` — blog posts
- `_projects/` — software/research projects (ScanNeo, ScanITD, ScanExitron, ...)
- `_bibliography/papers.bib` — publications
- `_pages/` — static pages (about, CV, books, blog index, ...)
- `_data/` — site data (socials, citations, CV, coauthors)

## Deployment

Pushing to `master` triggers [`.github/workflows/deploy.yml`](.github/workflows/deploy.yml), which builds the site and publishes it to GitHub Pages.

## License

Based on [al-folio](https://github.com/alshedivat/al-folio), available under the [MIT License](LICENSE).
