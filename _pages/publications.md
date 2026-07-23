---
layout: page
permalink: /publications/
title: Publications
description: Peer-reviewed scientific publications, book chapters, and journal articles.
nav: true
nav_order: 2
---

<!-- _pages/publications.md -->

<style type="text/css">
  .pubs-page {
    --pub-card-bg: var(--global-card-bg-color);
    --pub-muted: var(--global-text-color-light);
    --pub-border: var(--global-divider-color);
    --pub-link-bg: var(--global-code-bg-color);
    --pub-link-hover-bg: rgba(0, 113, 227, 0.08);
    max-width: 920px;
    margin: 0 auto;
    padding-bottom: 2rem;
  }

  body {
    padding-bottom: 0;
  }

  footer.fixed-bottom {
    position: static;
    margin-top: 2rem;
  }

  .pubs-heading {
    margin: 0 0 1.5rem;
    color: var(--pub-muted);
    font-size: 0.78rem;
    font-weight: 700;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    border-bottom: 1px solid var(--pub-border);
    padding-bottom: 0.5rem;
  }

  .pub-card {
    display: flex;
    align-items: center;
    gap: 1.75rem;
    padding: 1.4rem;
    margin-bottom: 1rem;
    background: var(--pub-card-bg);
    border: 1px solid var(--pub-border);
    border-radius: 8px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  }

  .pub-thumb {
    flex: 0 0 120px;
    width: 120px;
    height: 80px;
    overflow: hidden;
    border-radius: 6px;
    border: 1px solid var(--pub-border);
  }

  .pub-thumb img {
    width: 100%;
    height: 100%;
    display: block;
    object-fit: cover;
    transition: transform 0.35s ease;
  }

  .pub-card:hover .pub-thumb img {
    transform: scale(1.035);
  }

  .pub-body {
    flex: 1 1 auto;
    min-width: 0;
  }

  .pub-title {
    display: inline;
    color: var(--global-text-color);
    font-size: 1.06rem;
    font-weight: 700;
    line-height: 1.45;
  }

  .pub-title:hover {
    color: var(--global-theme-color);
    text-decoration: none;
  }

  .pub-authors {
    margin-top: 0.35rem;
    font-size: 0.92rem;
    color: var(--global-text-color);
  }

  .pub-venue {
    margin-top: 0.2rem;
    font-size: 0.88rem;
    color: var(--pub-muted);
    font-style: italic;
  }

  .pub-links {
    margin-top: 0.6rem;
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem;
  }

  .pub-link {
    display: inline-flex;
    align-items: center;
    gap: 0.35rem;
    padding: 0.25rem 0.6rem;
    background: var(--pub-link-bg);
    border: 1px solid var(--pub-border);
    border-radius: 4px;
    font-size: 0.82rem;
    color: var(--global-text-color) !important;
    text-decoration: none !important;
  }

  .pub-link:hover {
    background: var(--pub-link-hover-bg);
    border-color: var(--global-theme-color);
  }

  @media (max-width: 700px) {
    .pub-card {
      flex-direction: column;
      align-items: stretch;
      gap: 1rem;
      padding: 1rem;
    }

    .pub-thumb {
      flex-basis: auto;
      width: 100%;
      height: 150px;
    }
  }
</style>

<div class="pubs-page">
  <h2 class="pubs-heading">Publications</h2>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Tumor-specific lncRNA IGF1R-AS1 trans-regulates chromatin interactions associated with oncogenic MYC signaling preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.nature.com/articles/s41467-026-70814-4" target="_blank" rel="noopener noreferrer">Tumor-specific lncRNA IGF1R-AS1 trans-regulates chromatin interactions associated with oncogenic MYC signaling</a>
      <div class="pub-authors">Yongyong Yang#, <strong>Ting-You Wang</strong>#, Joshua Fry#, Yingming Li#, Qingshu Meng#, Nathan E. Patchen, Abhirami Ramakrishnan, et al.</div>
      <div class="pub-venue">Nature Communications, 2026, 17(1). (#Co-first authors)</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.nature.com/articles/s41467-026-70814-4" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="A genomic language model for chimera artifact detection in Nanopore direct RNA sequencing preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.nature.com/articles/s41467-026-68571-5" target="_blank" rel="noopener noreferrer">A genomic language model for chimera artifact detection in Nanopore direct RNA sequencing</a>
      <div class="pub-authors">Yangyang Li#, <strong>Ting-You Wang</strong>#, Qingxiang Guo, Yanan Ren, Xiaotong Lu, Qi Cao, Rendong Yang.</div>
      <div class="pub-venue">Nature Communications, 2026, 17(1), 1864. (#Co-first authors)</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.nature.com/articles/s41467-026-68571-5" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Ketone drink enhances therapeutic efficacy in prostate cancer by targeting EZH2 preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.nature.com/articles/s41389-025-00567-0" target="_blank" rel="noopener noreferrer">Ketone drink enhances therapeutic efficacy in prostate cancer by targeting EZH2</a>
      <div class="pub-authors">Chaehyun Yum, Richard Schaefer, Rui Wang, <strong>Ting-You Wang</strong>, Xiaotong Lu, Qi Liu, et al.</div>
      <div class="pub-venue">Oncogenesis, 2025, 14(1), 24.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.nature.com/articles/s41389-025-00567-0" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Androgen receptor-regulated lncRNA PRCAT71 promotes AR signaling through the interaction with KHSRP in prostate cancer preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.science.org/doi/10.1126/sciadv.adk6989" target="_blank" rel="noopener noreferrer">Androgen receptor-regulated lncRNA PRCAT71 promotes AR signaling through the interaction with KHSRP in prostate cancer</a>
      <div class="pub-authors">Yongyong Yang, <strong>Ting-You Wang</strong>, Qianru Li, Jiawen Lu, Yanan Ren, Adam B. Weiner, Joshua Fry, et al.</div>
      <div class="pub-venue">Science Advances, 2025, 11(15).</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.science.org/doi/10.1126/sciadv.adk6989" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="EZH2 directly methylates PARP1 and regulates its activity in cancer preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.science.org/doi/10.1126/sciadv.adl2804" target="_blank" rel="noopener noreferrer">EZH2 directly methylates PARP1 and regulates its activity in cancer</a>
      <div class="pub-authors">Qingshu Meng, Jiangchuan Shen, Yanan Ren, Qi Liu, Rui Wang, Qiaqia Li, Weihua Jiang, Quan Wang, Yixiang Zhang, Jonathan Trinidad, Xiaotong Lu, <strong>Tingyou Wang</strong>, Yanqiang Li, et al.</div>
      <div class="pub-venue">Science Advances, 2024, 10(48).</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.science.org/doi/10.1126/sciadv.adl2804" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Comparative analyses of gene networks mediating cancer metastatic potentials across lineage types preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://academic.oup.com/bib/article/25/4/bbae357/7714851" target="_blank" rel="noopener noreferrer">Comparative analyses of gene networks mediating cancer metastatic potentials across lineage types</a>
      <div class="pub-authors">Sheng Wang, Emily K Stroup, <strong>Ting-You Wang</strong>, Rendong Yang, Zhe Ji.</div>
      <div class="pub-venue">Briefings in Bioinformatics, 2024, 25(4).</div>
      <div class="pub-links">
        <a class="pub-link" href="https://academic.oup.com/bib/article/25/4/bbae357/7714851" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Chapter 5: Detecting Medium and Large Insertions and Deletions with transIndel preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://link.springer.com/protocol/10.1007/978-1-0716-2293-3_5" target="_blank" rel="noopener noreferrer">Chapter 5: Detecting Medium and Large Insertions and Deletions with transIndel</a>
      <div class="pub-authors"><strong>Ting-You Wang</strong> and Rendong Yang.</div>
      <div class="pub-venue">Variant Calling - Methods and Protocols, Humana, New York, NY, 2022, vol 2493.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://link.springer.com/protocol/10.1007/978-1-0716-2293-3_5" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Opposing transcriptional programs of KLF5 and AR emerge during therapy for advanced prostate cancer preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.nature.com/articles/s41467-021-26612-1" target="_blank" rel="noopener noreferrer">Opposing transcriptional programs of KLF5 and AR emerge during therapy for advanced prostate cancer</a>
      <div class="pub-authors">Meixia Che, Aashi Chaturvedi, Sarah Munro, Samuel Pitzen, Alexander Ling, Weijie Zhang, Joshua Mentzer, Sheng-Yu Ku, Loredana Puca, Yanyun Zhu, Andries Bergman, Tesa Severson, Colleen Forster, Yuzhen Liu, Jacob Hildebrand, Mark Daniel, <strong>Ting-You Wang</strong> et al.</div>
      <div class="pub-venue">Nature Communications, 2021, 12, 6377.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.nature.com/articles/s41467-021-26612-1" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="A pan-cancer transcriptome analysis of exitron splicing identifies novel cancer driver genes and neoepitopes preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00223-9" target="_blank" rel="noopener noreferrer">A pan-cancer transcriptome analysis of exitron splicing identifies novel cancer driver genes and neoepitopes</a>
      <div class="pub-authors"><strong>Ting-You Wang</strong>, Qi Liu, Yanan Ren, Sk. Kayum Alam, Li Wang, Zhu Zhu, Luke H. Hoeppner, Scott M. Dehm, Qi Cao, Rendong Yang.</div>
      <div class="pub-venue">Molecular Cell, 2021, 81:2246-2260 e2212. (Highlighted in Cancer Immunology Research, June 1, 2021)</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00223-9" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Integrated protocol for exitron and exitron-derived neoantigen identification using human RNA-seq data with ScanExitron and ScanNeo preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.sciencedirect.com/science/article/pii/S2666166721004949" target="_blank" rel="noopener noreferrer">Integrated protocol for exitron and exitron-derived neoantigen identification using human RNA-seq data with ScanExitron and ScanNeo</a>
      <div class="pub-authors"><strong>Ting-You Wang</strong> and Rendong Yang.</div>
      <div class="pub-venue">STAR Protocols, 2021, 2, 100788.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.sciencedirect.com/science/article/pii/S2666166721004949" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="LncGSEA: a versatile tool to infer lncRNA associated pathways from large-scale cancer transcriptome sequencing data preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://link.springer.com/article/10.1186/s12864-021-07900-y" target="_blank" rel="noopener noreferrer">LncGSEA: a versatile tool to infer lncRNA associated pathways from large-scale cancer transcriptome sequencing data</a>
      <div class="pub-authors">Yanan Ren, <strong>Ting-You Wang</strong>, Leah Anderton, Qi Cao, Rendong Yang.</div>
      <div class="pub-venue">BMC Genomics, 2021, 22(1):574.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://link.springer.com/article/10.1186/s12864-021-07900-y" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Identification of 38 novel loci for systemic lupus erythematosus and genetic heterogeneity between ancestral groups preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.nature.com/articles/s41467-021-21049-y" target="_blank" rel="noopener noreferrer">Identification of 38 novel loci for systemic lupus erythematosus and genetic heterogeneity between ancestral groups</a>
      <div class="pub-authors">Yong-Fei Wang, Yan Zhang, Zhiming Lin, Huoru Zhang, <strong>Ting-You Wang</strong>, Yujie Cao, David L Morris, et al.</div>
      <div class="pub-venue">Nature Communications, 2021, 12, 772.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.nature.com/articles/s41467-021-21049-y" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://academic.oup.com/gigascience/article/9/8/giaa089/5898622" target="_blank" rel="noopener noreferrer">ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation</a>
      <div class="pub-authors"><strong>Ting-You Wang</strong>, Rendong Yang.</div>
      <div class="pub-venue">GigaScience, 2020, 9.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://academic.oup.com/gigascience/article/9/8/giaa089/5898622" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="RAS internal tandem duplication disrupts GTPase-activating protein (GAP) binding to activate oncogenic signaling preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.jbc.org/article/S0021-9258(17)48957-5/fulltext" target="_blank" rel="noopener noreferrer">RAS internal tandem duplication disrupts GTPase-activating protein (GAP) binding to activate oncogenic signaling</a>
      <div class="pub-authors">Andrew C Nelson, Thomas J Turbyville, Srisathiyanarayanan Dharmaiah, Megan Rigby, Rendong Yang, <strong>Ting-You Wang</strong>, John Columbus, et al.</div>
      <div class="pub-venue">Journal of Biological Chemistry, 2020, 295, 9335-9348.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.jbc.org/article/S0021-9258(17)48957-5/fulltext" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Identification of regulatory modules that stratify lupus disease mechanism through integrating multi-omics data preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.sciencedirect.com/science/article/pii/S2162253119303750" target="_blank" rel="noopener noreferrer">Identification of regulatory modules that stratify lupus disease mechanism through integrating multi-omics data</a>
      <div class="pub-authors"><strong>Ting-You Wang</strong>, Yong-Fei Wang, Yan Zhang, Jiangshan Jane Shen, Mengbiao Guo, Jing Yang, et al.</div>
      <div class="pub-venue">Molecular Therapy - Nucleic Acids, 2020, 19, 318-329.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.sciencedirect.com/science/article/pii/S2162253119303750" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="ScanNeo: identifying indel derived neoantigens using RNA-Seq data preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://academic.oup.com/bioinformatics/article/35/20/4159/5382215" target="_blank" rel="noopener noreferrer">ScanNeo: identifying indel derived neoantigens using RNA-Seq data</a>
      <div class="pub-authors"><strong>Ting-You Wang</strong>, Li Wang, Sk. Kayum Alam, Luke H. Hoeppner, Rendong Yang.</div>
      <div class="pub-venue">Bioinformatics, 2019, 35(20), 4159-4161.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://academic.oup.com/bioinformatics/article/35/20/4159/5382215" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="HLA-IMPUTER: an easy to use web application for HLA imputation and association analysis using population-specific reference panels preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://academic.oup.com/bioinformatics/article/35/7/1244/5101258" target="_blank" rel="noopener noreferrer">HLA-IMPUTER: an easy to use web application for HLA imputation and association analysis using population-specific reference panels</a>
      <div class="pub-authors">Jiangshan Shen, Chao Yang, Yong-Fei Wang, <strong>Ting-You Wang</strong>, Mengbiao Guo, et al.</div>
      <div class="pub-venue">Bioinformatics, 2018, 35(7), 1244-1246.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://academic.oup.com/bioinformatics/article/35/7/1244/5101258" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Identification of ST3AGL4, MFHAS1, CSNK2A2 and CD226 as loci associated with systemic lupus erythematosus (SLE) and evaluation of SLE genetics in drug repositioning preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://ard.bmj.com/content/77/7/1078" target="_blank" rel="noopener noreferrer">Identification of ST3AGL4, MFHAS1, CSNK2A2 and CD226 as loci associated with systemic lupus erythematosus (SLE) and evaluation of SLE genetics in drug repositioning</a>
      <div class="pub-authors">Yong-Fei Wang, Yan Zhang, Zhengwei Zhu, <strong>Ting-You Wang</strong>, David L Morris, Jiangshan Jane Shen, et al.</div>
      <div class="pub-venue">Annals of the Rheumatic Diseases, 2018, 77(7):1078-1084.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://ard.bmj.com/content/77/7/1078" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Regulatory and evolutionary signatures of sex-biased genes on both the X chromosome and the autosomes preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://link.springer.com/article/10.1186/s13293-017-0156-4" target="_blank" rel="noopener noreferrer">Regulatory and evolutionary signatures of sex-biased genes on both the X chromosome and the autosomes</a>
      <div class="pub-authors">Jiangshan Shen, <strong>Ting-You Wang</strong>, Wanling Yang.</div>
      <div class="pub-venue">Biology of Sex Differences, 2017, 8(1):35.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://link.springer.com/article/10.1186/s13293-017-0156-4" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Genome-wide association meta-analysis in Chinese and European individuals identifies ten new loci associated with systemic lupus erythematosus preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://www.nature.com/articles/ng.3603" target="_blank" rel="noopener noreferrer">Genome-wide association meta-analysis in Chinese and European individuals identifies ten new loci associated with systemic lupus erythematosus</a>
      <div class="pub-authors">David L Morris, Yujun Sheng, Yan Zhang, Yong-Fei Wang, Zhengwei Zhu, Philip Tombleson, Lingyan Chen, Deborah S Cunninghame Graham, James Bentham, Amy L Roberts, Ruoyan Chen, Xianbo Zuo, <strong>Tingyou Wang</strong>, Leilei Wen, Chao Yang, Lu Liu, Lulu Yang, Feng Li, Yuanbo Huang, Xianyong Yin, Sen Yang, Lars Rönnblom, Barbara G Fürnrohr, Reinhard E Voll, Georg Schett et al.</div>
      <div class="pub-venue">Nature Genetics, 2016, 48, 940-946.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://www.nature.com/articles/ng.3603" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Meta-analysis of GWAS on two Chinese populations followed by replication identifies novel genetic variants on the X chromosome associated with systemic lupus erythematosus preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://academic.oup.com/hmg/article/24/1/274/2900742" target="_blank" rel="noopener noreferrer">Meta-analysis of GWAS on two Chinese populations followed by replication identifies novel genetic variants on the X chromosome associated with systemic lupus erythematosus</a>
      <div class="pub-authors">Yan Zhang, Jing Zhang, Jing Yang, Yongfei Wang, Lu Zhang, Xianbao Zuo, Liangdan Sun, Hai-Feng Pan, Nattiya Hirankarn, <strong>Tingyou Wang</strong>, Ruoyan Chen, Dingge Ying, et al.</div>
      <div class="pub-venue">Human Molecular Genetics, 2015, 24(1), 274-284.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://academic.oup.com/hmg/article/24/1/274/2900742" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Chromosome-8-Coded Proteome of Chinese Chromosome Proteome Data Set (CCPD) 2.0 with Partial Immunohistochemical Verifications preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://pubs.acs.org/doi/10.1021/pr400902u" target="_blank" rel="noopener noreferrer">Chromosome-8-Coded Proteome of Chinese Chromosome Proteome Data Set (CCPD) 2.0 with Partial Immunohistochemical Verifications</a>
      <div class="pub-authors">Yang Liu, Wantao Ying, Zhe Ren, Wei Gu, Yang Zhang, Guoquan Yan, Pengyuan Yang, Yinkun Liu, Xuefei Yin, Cheng Chang, Jing Jiang, Fengxu Fan, Chengpu Zhang, Ping Xu, Quanhui Wang, Bo Wen, Liang Lin, <strong>Tingyou Wang</strong>, Chaoqin Du et al.</div>
      <div class="pub-venue">Journal of Proteome Research, 2014, 13(1), 126-136.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://pubs.acs.org/doi/10.1021/pr400902u" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="Qualitative and quantitative expression status of the human chromosome 20 genes in cancer tissues and the representative cell lines preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://pubs.acs.org/doi/10.1021/pr3007909" target="_blank" rel="noopener noreferrer">Qualitative and quantitative expression status of the human chromosome 20 genes in cancer tissues and the representative cell lines</a>
      <div class="pub-authors">Quanhui Wang#, Bo Wen#, Guangrong Yan#, Junying Wei#, Liqi Xie#, Shaohang Xu, Dahai Jiang, <strong>Tingyou Wang</strong>, Ruo Zhou, Haiyi Zhao et al.</div>
      <div class="pub-venue">Journal of Proteome Research, 2013, 12(1), 151-161. (#Co-first authors)</div>
      <div class="pub-links">
        <a class="pub-link" href="https://pubs.acs.org/doi/10.1021/pr3007909" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>

  <article class="pub-card">
    <div class="pub-thumb">
      <img src="/assets/img/publication_preview/default_paper.jpg" alt="A predicted protein-protein interaction network of the filamentous fungus Neurospora crassa preview" />
    </div>
    <div class="pub-body">
      <a class="pub-title" href="https://pubs.rsc.org/en/content/articlelanding/2011/mb/c1mb05028a" target="_blank" rel="noopener noreferrer">A predicted protein-protein interaction network of the filamentous fungus Neurospora crassa</a>
      <div class="pub-authors"><strong>Ting-You Wang</strong>, Fei He, Qi-Wen Hu, Ziding Zhang.</div>
      <div class="pub-venue">Molecular BioSystems, 2011, 7, 2278-2285.</div>
      <div class="pub-links">
        <a class="pub-link" href="https://pubs.rsc.org/en/content/articlelanding/2011/mb/c1mb05028a" target="_blank" rel="noopener noreferrer"><i class="fa-solid fa-file-pdf"></i>Paper</a>
      </div>
    </div>
  </article>
</div>