---
layout: post
title: "Rethinking Bioinformatics Expertise in the Era of AI: A Newbie's Journey as a Custodian"
date: 2026-07-23
category: opinion
comments: true
tags: [Bioinformatics, Artificial Intelligence, Agentic Coding, Antigravity, Personal Perspective]
image: /figures/2026-07-23-rethinking-bioinformatics-expertise/cover.jpg
---

![Bioinformatics AI Synergy](/figures/2026-07-23-rethinking-bioinformatics-expertise/cover.jpg)

Earlier this year, in 2026, I started bringing advanced AI agents into my daily coding and research workflows. Tools like Claude Code and Google Antigravity became part of my regular toolkit. If I'm being honest, I am still very much a newbie at using them—I am still wrapping my head around advanced features like Antigravity **SKILLS** and modular agent configurations. 

Yet, even as a beginner, this transition has forced me to grapple with a fundamental question: Is AI a threat to our careers as bioinformaticians, or is it the ultimate tool for augmenting our abilities?

The answer is both. But the key to surviving—and thriving—lies in one word: **testing**.

<!--more-->

My personal journey over the last few months matches up perfectly with a recent perspective in *npj Digital Medicine* titled ["Rethinking bioinformatics expertise in the era of artificial intelligence"](https://www.nature.com/articles/s41746-026-02777-1) by Goh et al. (May 2026). Reading their analysis felt like reading a roadmap for exactly what I’ve been experiencing at my desk. 

Here is my personal perspective on why the bioinformatics profession isn't dying, but why we must urgently transition from being passive consumers of AI to active custodians who know how to test, validate, and guide these systems.

---

## 1. The Real Threat: The "Appearance" of Science

When I first booted up Claude Code and Antigravity, I felt a familiar pang of AI anxiety. The threat felt real: when an LLM can write an entire differential expression pipeline or generate a structural model in seconds, it’s easy to feel obsolete. 

But as I began using these tools for actual work, the nature of that threat shifted. The danger isn’t that AI is going to replace us because it’s "smarter" than us. The danger is that AI is incredibly good at democratizing the **appearance of science**. 

An agent can confidently write a script, generate a beautiful plot, and output a highly plausible-sounding biological interpretation. But the agent doesn't actually understand biology. It doesn't know the physical reality of a sequencer, a batch effect, or an intrinsically disordered protein region. It is only matching statistical patterns in its training data. 

If a non-expert runs the model, they will get a clean, professional-looking report. But they will have no way of knowing if the analysis is biologically nonsensical or confidently wrong. The threat of AI isn't automation; it's the proliferation of high-tech noise masquerading as scientific truth.

---

## 2. What Daily Work Taught Me: Testing is Everything

If there is one thing my rookie journey with Google Antigravity and Claude Code has taught me, it is this: **testing is the most critical skill we have in the AI era.**

When you delegate code generation or pipeline execution to an AI agent, you aren't offloading the responsibility of the science. If anything, your responsibility multiplies. You are now the gatekeeper who must verify that the code behaves correctly under the unique constraints of your experiment.

In my daily work, I’ve realized that I spend less time typing boilerplate code and significantly more time designing test cases:
* Is the generated script using the correct reference genome build (e.g., GRCh38)?
* Does the population stratification correction in this automated GWAS pipeline actually work on our specific cohort?
* Are the library prep strand-specificity parameters correctly specified in the RNA-seq commands?

If you don't know how to write robust tests, if you don't understand the underlying biology well enough to recognize a subtle statistical artifact, then AI is a threat. But if you have the domain expertise to design rigorous validation checkpoints, AI becomes an extraordinary accelerant. It handles the manual drafting, leaving you free to focus on experimental design and validation.

```
[AI generates draft pipeline] ────> [Bioinformatician designs biological tests] ────> [Validated Discovery]
                                                      │
                                           (Catches silent bugs,
                                          biases, & hallucinations)
```

---

## 3. Transitioning from Consumer to Custodian

In their paper, Goh et al. talk about shifting from an "AI Consumer" to an "AI Custodian." This concept clicked for me as I explored Antigravity's **SKILLS** framework. 

As a newbie, it’s tempting to just write a simple prompt and copy-paste the output. But being a custodian means treating AI prompts and workflows as living methodological assets. It means:
1. **Structuring Prompts and Agents:** Crafting domain-specific, expert-curated instructions (such as Chain-of-Thought prompts) that steer LLMs toward biological reality.
2. **Auditing the Data:** Realizing that model performance depends entirely on data quality. We must identify batch effects, label noise, and metadata standards before the data ever touches a model.
3. **Interrogating Explanations:** Refusing to accept explainability maps (like SHAP values) at face value. We must test whether the features a classifier relies on are biologically mechanistic or just housekeeping genes correlated with technical noise.

```
                    AI CONSUMER                             AI CUSTODIAN
             ┌─────────────────────────┐             ┌─────────────────────────┐
             │  • Blindly trusts code  │             │  • Writes rigorous tests│
             │  • Copies raw outputs   │     VS.     │  • Audits training data │
             │  • Ignores edge cases   │             │  • Evaluates XAI maps   │
             └─────────────────────────┘             └─────────────────────────┘
```

---

## 4. Bridging the Vocabulary Gap in Digital Medicine

Another area where our expertise is irreplaceable is in managing what the authors call the **vocabulary gap**. 

In my projects, I've noticed a persistent communication disconnect:
* **The Computational Side** talks about *pipelines, database schemas, features, and model weights*.
* **The Clinical/Life Science Side** talks about *experimental assays, clinical provenance, biomarker validation, and patient outcomes*.

A machine learning engineer might see a clinical dataset as a matrix of numbers to optimize. But a bioinformatician knows the messy reality behind those numbers—how the samples were collected, which assays are sensitive to temperature, and how a change in hospital EHR software can introduce silent distribution shifts. 

As bioinformaticians, we act as the translators. We understand enough of the computer science to configure agents and build code, and enough of the biology to know if the clinical assumptions make sense. No generic AI tool can bridge that gap; it requires hard-won, interdisciplinary human experience.

---

## 5. The Institutional Challenge: Ethics and Regulations

As I look ahead, I realize that the challenge isn't just about polishing our individual scripts. It's about how our institutions deploy AI. 

We are entering an era of heavy regulation. The **EU AI Act**, which went into effect in August 2024, classifies AI clinical decision support systems as **high-risk**. By **August 2027**, institutions will be legally required to provide data governance, human oversight, subgroup performance testing, and continuous post-market monitoring.

This means we must treat ethics and bias auditing not as compliance paperwork, but as core scientific rigor. 

For instance, if a diagnostic model is trained primarily on one demographic, it will fail on others. If we deploy it blindly, we aren't just doing bad science; we are actively endangering patients. Bioinformaticians are the ones who must write the subgroup audits, stratify model performance by ancestry or age, and flag where models underperform.

> **My Personal "In Practice" Rule**
> Every time I use an AI agent to build a translational model, I run a subgroup stratification audit. If the model drops in accuracy for a specific demographic, I document it and flag it as a limitation. It’s not extra work; it is the work.

---

## Conclusion: The Road Ahead

Six months ago, I was worried that developer agents would make my coding skills obsolete. Today, I see things differently. 

Yes, AI is a threat if we remain passive consumers who copy and paste code without understanding it. But if we step up as custodians—if we focus on designing tests, curating clean biological data, translating across disciplines, and auditing for bias—AI becomes the most powerful augmentant we have ever seen.

I am still a newbie, and I still have a lot to learn about optimizing SKILLS and orchestrating agents. But the most important lesson I’ve learned doesn't require being an AI wizard. It just requires being a good scientist: **Never trust the machine blindly. Design the test, verify the biology, and be the expert behind the keyboard.**

***

*What are your experiences using developer agents like Claude Code or Antigravity in your research? Let's discuss in the comments below! For the original perspective, read the full paper in npj Digital Medicine: [Rethinking bioinformatics expertise in the era of artificial intelligence](https://www.nature.com/articles/s41746-026-02777-1).*
