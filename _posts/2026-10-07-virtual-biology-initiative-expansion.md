---
layout: post
title: "$1.8 Billion for AI-Ready Biology: A Bioinformatician's Take on the Virtual Biology Initiative Expansion"
date: 2026-10-07
category: opinion
comments: true
tags: [Virtual Cell, Artificial Intelligence, Single-Cell, Open Data, Biohub, Personal Perspective]
image: /figures/2026-10-07-virtual-biology-initiative-expansion/cover.jpg
---

![Virtual Biology Initiative cover](/figures/2026-10-07-virtual-biology-initiative-expansion/cover.jpg)

Today Biohub, the U.S. Department of Energy, the National Institutes of Health, and a group of industry partners announced a **\\$1.8 billion** commitment to generate and openly share biological data built for training predictive AI models. They call it the largest coordinated investment in AI-ready biological data to date.

My first reaction was a mix of excitement and a familiar question: is the bottleneck for AI in biology really models, or has it been the data all along?

This announcement is a strong bet on the second answer.

<!--more-->

The full release is on the Biohub site: ["AI-ready biological data: \\$1.8 billion global commitment"](https://biohub.org/news/virtual-biology-initiative-expansion/). It expands the [Virtual Biology Initiative](https://biohub.org/news/virtual-biology-initiative/), first announced in April 2026, whose long-term goal is a "virtual cell": a model that can predict how cells respond to perturbations, drugs, and disease.

Below is my personal read on what is in the package, why it matters, and what I'll be watching for as someone who works with this kind of data every day.

---

## 1. What Was Actually Committed

The headline number is \\$1.8 billion, but it is made of quite different kinds of contributions:

- **Biohub — \\$500 million (founding commitment).** About \\$400M goes to new measurement technology: cryo-electron tomography at near-atomic resolution inside cells, microscopy that images millions to billions of cells in living tissue, and engineering tools to build and perturb biology from molecules up to whole organisms. About \\$100M funds research outside Biohub.
- **U.S. Department of Energy — over \\$500 million over five years**, delivered through DOE's cross-agency **Genesis Mission**. This brings exascale supercomputing, X-ray and neutron scattering, cryo-EM, autonomous labs, and facilities such as the Joint Genome Institute and the Environmental Molecular Sciences Laboratory.
- **NIH — coordination of over \\$500 million in prior federal investment** via the [Bio Genesis Mission](https://www.nih.gov/bio-genesismission). Rather than new money, this is about making existing datasets, repositories (NCBI, NLM-catalogued resources) and Common Fund atlases interoperable and AI-ready.
- **Google DeepMind, Isomorphic Labs, and Meta — \\$300 million combined** for multimodal datasets and technologies.
- **NVIDIA** contributes accelerated compute, domain software and expertise; **Renaissance Philanthropy** helps grow data-generation funding.

The scientific partners include the Allen Institute, Broad Institute, Gladstone Institutes, Human Cell Atlas, Human Protein Atlas, and Wellcome Sanger Institute.

---

## 2. Data Is the Real Bottleneck

Over the past two years, we've seen a wave of single-cell "foundation models." Many are impressive on paper. But when you benchmark them carefully against simple baselines for perturbation prediction, the gains often shrink or vanish.

In my view, a big reason is the training data. Most public single-cell data are observational atlases: they tell you what cells _look like_, not how they _respond_ when you push on them. Perturbation data exist, but they cover a narrow slice of cell types, conditions, and readouts, and they come from labs using different protocols.

That is why the stated goals of this expansion caught my attention:

1. **Expand cell response data** across far more cell types and conditions than have been studied.
2. **Build and validate measurement technologies** at greater scale, speed, and accuracy.
3. **Create shared standards, common identifiers, and a single access point**, so datasets from different partners actually work together.

Point 3 is the least glamorous and, I suspect, the most important.

---

## 3. Standards Are Where the Hard Work Lives

Anyone who has tried to merge two public scRNA-seq datasets knows the pain: different gene identifiers, different reference annotations, missing metadata on donor, tissue, dissociation protocol, or chemistry version, and batch effects that swamp biology.

Biohub has a real track record here. [CELLxGENE](https://cellxgene.cziscience.com/) standardized schemas and ontology terms across thousands of datasets, and the [CryoET Data Portal](https://cryoetdataportal.czscience.com/) did similar work for tomography. Projects like [Tabula Sapiens](https://tabula-sapiens.sf.czbiohub.org/), [OpenCell](https://opencell.sf.czbiohub.org/), and [Zebrahub](https://zebrahub.sf.czbiohub.org/) showed what consistent, well-annotated reference data can look like.

Scaling that discipline across government labs, NIH repositories, and industry partners is a much bigger challenge. A "single access point" is only as useful as the metadata behind it. If this initiative gets the schemas, controlled vocabularies, and QC standards right, it could matter for years, regardless of which model architecture wins.

> **My Personal "In Practice" Rule**
> Before I train or fine-tune anything on public data, I spend time on the metadata: checking gene ID versions, harmonizing cell-type labels to an ontology, and looking for batch structure. AI-ready data means this work is done _before_ the data are released, not left to every downstream user.

---

## 4. Open Data, Industry Money: Questions Worth Asking

The release emphasizes open sharing. Pushmeet Kohli of Google DeepMind put it plainly: _"We will not solve this challenge without open, experimental biological data at an unprecedented scale"_.

I agree, and I'm glad to see industry money going into public data generation. Still, there are questions I'd like answered as details emerge:

- **Access terms.** Will the data be fully open (like CELLxGENE), or under controlled access? Human data will reasonably need protections; how will that balance be struck?
- **Timing.** Will industry partners get early access before public release?
- **Benchmarks.** Will there be held-out, community-run benchmarks so that "virtual cell" claims can be tested fairly, rather than each group reporting on its own splits?
- **Diversity.** Will cell lines and donor samples reflect human genetic diversity, or repeat the ancestry bias already present in many genomics resources?

None of these are reasons for cynicism. They are the details that decide whether a resource becomes a true public good.

---

## 5. What It Means for the Rest of Us

Most of us won't be generating petabytes of perturbation data or running exascale jobs. But we will be the users and, hopefully, the critics of what comes out.

For working bioinformaticians, I see three practical takeaways:

1. **Learn the standards early.** Get comfortable with AnnData/CELLxGENE schemas, cell ontologies, and the tooling around them (e.g., `cellxgene-census`). Data that follow these conventions will be the easiest to combine with this new resource.
2. **Benchmark before you believe.** When "virtual cell" models trained on this data appear, test them on your own systems and against simple baselines.
3. **Contribute metadata, not just data.** Well-annotated smaller datasets will be more valuable to this ecosystem than large poorly annotated ones.

---

## Conclusion: A Bet on Infrastructure

As Alex Rives, Biohub's Head of Science, said, _"the creation of a virtual cell is one of the most important challenges for the next era of science."_ I share the excitement, but what makes me optimistic about this announcement is not the virtual cell goal itself. It's that so much of the money is going into measurement technology, standards, and shared infrastructure rather than another round of model building on the same limited data.

If the initiative delivers truly open, well-annotated, interoperable perturbation data at scale, it will help everyone, including people who are skeptical of today's foundation models. **Better data first, then better models.**

---

_What do you think: is data or modeling the main bottleneck for virtual cell models? Let's discuss in the comments below! Read the original announcement on the Biohub website: [AI-ready biological data: \\$1.8 billion global commitment](https://biohub.org/news/virtual-biology-initiative-expansion/)._
