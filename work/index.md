---
title: Work
layout: page
---

<style>
  /* Premium Work Page Stylesheet */
  .work-container {
    color: #2d3748;
    font-family: 'Open Sans', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    line-height: 1.6;
  }

  /* Header Overrides & Styling */
  .work-header {
    text-align: center;
    margin-bottom: 2.5rem;
    padding-bottom: 1.5rem;
    border-bottom: 1px solid #edf2f7;
  }
  
  .work-title {
    font-family: 'Vollkorn', Georgia, serif;
    font-size: 2.4rem;
    color: #1a202c;
    margin-bottom: 0.5rem;
    font-weight: 500;
  }

  .work-subtitle {
    font-size: 1.1rem;
    color: #718096;
    max-width: 600px;
    margin: 0 auto;
    font-weight: 300;
  }

  .work-actions {
    display: flex;
    justify-content: center;
    gap: 12px;
    margin-top: 1.2rem;
  }

  /* Profile Buttons */
  .profile-btn {
    display: inline-flex;
    align-items: center;
    gap: 8px;
    padding: 8px 16px;
    border-radius: 8px;
    font-size: 0.9rem;
    font-weight: 600;
    text-decoration: none !important;
    transition: all 0.2s ease;
    border: 1px solid #e2e8f0;
    background-color: #fff;
    color: #4a5568 !important;
  }

  .profile-btn:hover {
    transform: translateY(-2px);
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
    border-color: #3aafc1;
    color: #3aafc1 !important;
  }

  .profile-btn svg {
    width: 16px;
    height: 16px;
    fill: currentColor;
  }

  /* Stat Dashboard Grid */
  .stats-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 15px;
    margin-bottom: 3rem;
  }

  .stat-card {
    background: linear-gradient(135deg, #ffffff 0%, #f7fafc 100%);
    border: 1px solid #e2e8f0;
    border-radius: 12px;
    padding: 20px 15px;
    text-align: center;
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.02);
    transition: all 0.3s ease;
  }

  .stat-card:hover {
    transform: translateY(-3px);
    box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.05);
    border-color: #3aafc1;
  }

  .stat-num {
    font-size: 2.2rem;
    font-weight: 700;
    color: #3aafc1;
    line-height: 1;
    margin-bottom: 5px;
  }

  .stat-label {
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: #718096;
    font-weight: 600;
  }

  /* Section Titles */
  .section-title {
    font-family: 'Vollkorn', Georgia, serif;
    font-size: 1.65rem;
    color: #1a202c;
    margin-top: 2.5rem;
    margin-bottom: 1.2rem;
    display: flex;
    align-items: center;
    gap: 10px;
  }

  .section-title::after {
    content: '';
    flex-grow: 1;
    height: 1px;
    background-color: #e2e8f0;
  }

  .section-title-icon {
    color: #3aafc1;
    display: flex;
    align-items: center;
  }

  /* Research Areas */
  .research-grid {
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 20px;
    margin-bottom: 2.5rem;
  }

  .research-card {
    background-color: #fff;
    border: 1px solid #edf2f7;
    border-radius: 12px;
    padding: 20px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.01);
    transition: all 0.2s ease;
  }

  .research-card:hover {
    border-color: #3aafc1;
    box-shadow: 0 6px 12px rgba(58, 175, 193, 0.05);
  }

  .research-title {
    font-size: 1.1rem;
    font-weight: 600;
    color: #2d3748;
    margin-bottom: 10px;
    display: flex;
    align-items: center;
    gap: 8px;
  }

  .research-title svg {
    color: #3aafc1;
    width: 18px;
    height: 18px;
  }

  .research-desc {
    font-size: 0.9rem;
    color: #4a5568;
    margin: 0;
  }

  /* Software Grid & Cards */
  .software-grid {
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 20px;
    margin-bottom: 3rem;
  }

  .software-card {
    background-color: #fff;
    border: 1px solid #e2e8f0;
    border-radius: 12px;
    padding: 22px;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.01);
    transition: all 0.3s cubic-bezier(0.25, 0.8, 0.25, 1);
    display: flex;
    flex-direction: column;
    justify-content: space-between;
  }

  .software-card:hover {
    transform: translateY(-5px);
    border-color: #3aafc1;
    box-shadow: 0 12px 20px rgba(58, 175, 193, 0.08);
  }

  .card-top {
    margin-bottom: 15px;
  }

  .card-header {
    display: flex;
    justify-content: space-between;
    align-items: flex-start;
    margin-bottom: 8px;
  }

  .card-title {
    font-size: 1.25rem;
    font-weight: 600;
    margin: 0;
  }

  .card-title a {
    color: #1a202c !important;
    text-decoration: none !important;
    transition: color 0.2s;
  }

  .card-title a:hover {
    color: #3aafc1 !important;
  }

  .card-github-icon {
    color: #718096;
    transition: color 0.2s;
  }

  .card-github-icon:hover {
    color: #1a202c;
  }

  .card-desc {
    font-size: 0.9rem;
    color: #4a5568;
    margin: 0 0 12px 0;
  }

  .tag-container {
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
  }

  .tag {
    font-size: 0.72rem;
    padding: 2px 8px;
    border-radius: 9999px;
    font-weight: 500;
    background-color: #edf2f7;
    color: #4a5568;
  }

  .tag.primary {
    background-color: rgba(58, 175, 193, 0.1);
    color: #0a9cb2;
  }

  .card-action-link {
    font-size: 0.85rem;
    font-weight: 600;
    color: #3aafc1 !important;
    text-decoration: none !important;
    display: inline-flex;
    align-items: center;
    gap: 4px;
    transition: color 0.2s;
  }

  .card-action-link:hover {
    color: #0a9cb2 !important;
  }

  /* Publications List */
  .pub-list {
    display: flex;
    flex-direction: column;
    gap: 15px;
    margin-bottom: 3rem;
  }

  .pub-item {
    background-color: #ffffff;
    border-left: 4px solid #3aafc1;
    padding: 15px 20px;
    border-radius: 0 10px 10px 0;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.02);
    border-top: 1px solid #edf2f7;
    border-right: 1px solid #edf2f7;
    border-bottom: 1px solid #edf2f7;
    transition: all 0.2s ease;
  }

  .pub-item:hover {
    background-color: #f7fafc;
    border-left-color: #0a9cb2;
    transform: translateX(4px);
  }

  .pub-title {
    font-weight: 600;
    font-size: 1rem;
    color: #1a202c;
    margin-bottom: 6px;
    line-height: 1.4;
  }

  .pub-title a {
    color: #1a202c !important;
    text-decoration: none !important;
  }

  .pub-title a:hover {
    color: #3aafc1 !important;
  }

  .pub-authors {
    font-size: 0.88rem;
    color: #4a5568;
    margin-bottom: 6px;
  }

  .pub-authors strong {
    color: #1a202c;
    font-weight: 600;
  }

  .pub-footer {
    display: flex;
    justify-content: space-between;
    align-items: center;
    font-size: 0.82rem;
    flex-wrap: wrap;
    gap: 10px;
  }

  .pub-journal {
    color: #718096;
    font-style: italic;
  }

  .pub-stats {
    display: flex;
    align-items: center;
    gap: 10px;
  }

  .pub-citation-badge {
    background-color: rgba(58, 175, 193, 0.1);
    color: #0a9cb2;
    padding: 2px 8px;
    border-radius: 9999px;
    font-weight: 600;
    font-size: 0.78rem;
  }

  .pub-link {
    color: #3aafc1 !important;
    text-decoration: none !important;
    font-weight: 600;
    display: inline-flex;
    align-items: center;
    gap: 3px;
  }

  .pub-link:hover {
    color: #0a9cb2 !important;
  }

  /* Responsive Adjustments */
  @media (max-width: 768px) {
    .stats-grid {
      grid-template-columns: repeat(2, 1fr);
    }
    .research-grid {
      grid-template-columns: 1fr;
    }
    .software-grid {
      grid-template-columns: 1fr;
    }
  }

  @media (max-width: 480px) {
    .stats-grid {
      grid-template-columns: 1fr;
    }
    .work-title {
      font-size: 2rem;
    }
  }
</style>

<div class="work-container">
  
  <!-- Page Header -->
  <div class="work-header">
    <h2 class="work-title">Research & Code</h2>
    <p class="work-subtitle">
      Bioinformatics Research Associate at Northwestern University. Deciphering cancer regulatory complexity, lupus genetics, and transcriptomics through robust open-source systems.
    </p>
    <div class="work-actions">
      <a href="https://scholar.google.com/citations?user=FVx_uDIAAAAJ&hl=en" class="profile-btn" target="_blank" rel="noopener noreferrer">
        <svg viewBox="0 0 24 24"><path d="M12 2L1 9l11 7 9-5.73V17h2V9L12 2zm5.12 9.24L12 14.5l-5.12-3.26L12 8.35l5.12 2.89zM12 17.5l-4.5-2.86v3.25L12 20.8l4.5-2.91v-3.25l-4.5 2.86z"/></svg>
        Google Scholar
      </a>
      <a href="https://github.com/dolittle007" class="profile-btn" target="_blank" rel="noopener noreferrer">
        <svg viewBox="0 0 24 24"><path d="M12 2A10 10 0 0 0 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.9-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.9 1.52 2.34 1.07 2.91.83.1-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0 0 12 2z"/></svg>
        GitHub Profile
      </a>
    </div>
  </div>

  <!-- Stats Dashboard -->
  <div class="stats-grid">
    <div class="stat-card">
      <div class="stat-num">1,155+</div>
      <div class="stat-label">Total Citations</div>
    </div>
    <div class="stat-card">
      <div class="stat-num">14</div>
      <div class="stat-label">h-index</div>
    </div>
    <div class="stat-card">
      <div class="stat-num">16</div>
      <div class="stat-label">i10-index</div>
    </div>
    <div class="stat-card">
      <div class="stat-num">21</div>
      <div class="stat-label">Publications</div>
    </div>
  </div>

  <!-- Research Themes -->
  <div class="section-title">
    <span class="section-title-icon">
      <svg width="24" height="24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-compass"><circle cx="12" cy="12" r="10"></circle><polygon points="16.24 7.76 14.12 14.12 7.76 16.24 9.88 9.88 16.24 7.76"></polygon></svg>
    </span>
    Research Focus Areas
  </div>
  
  <div class="research-grid">
    <div class="research-card">
      <div class="research-title">
        <svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M4.5 16.5c-1.5 1.26-2.5 3.19-2.5 5.5h20c0-2.31-1-4.24-2.5-5.5M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8z"/></svg>
        Lupus Susceptibility & Multi-Omics
      </div>
      <p class="research-desc">
        Investigating Systemic Lupus Erythematosus (SLE) etiology via large-scale GWAS meta-analyses across diverse cohorts, resolving ancestral heterogeneity, and mapping molecular risk cascades.
      </p>
    </div>
    <div class="research-card">
      <div class="research-title">
        <svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="12 2 2 7 12 12 22 7 12 2"></polygon><polyline points="2 17 12 22 22 17"></polyline><polyline points="2 12 12 17 22 12"></polyline></svg>
        Alternative Splicing & Splicing Drivers
      </div>
      <p class="research-desc">
        Developing computational methods to capture aberrant non-canonical splicing events, such as exitron splicing, and examining their functional roles as neoepitopes and drivers in pan-cancer cohorts.
      </p>
    </div>
    <div class="research-card">
      <div class="research-title">
        <svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><rect x="2" y="2" width="20" height="8" rx="2" ry="2"></rect><rect x="2" y="14" width="20" height="8" rx="2" ry="2"></rect><line x1="6" y1="6" x2="6.01" y2="6"></line><line x1="6" y1="18" x2="6.01" y2="18"></line></svg>
        Immunogenomics & Neoantigen Prediction
      </div>
      <p class="research-desc">
        Constructing algorithmic pipelines (e.g., ScanNeo) to robustly predict frame-shifting, indel-derived neoantigens from transcriptome sequencing, aiding personalized cancer immunotherapies.
      </p>
    </div>
    <div class="research-card">
      <div class="research-title">
        <svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5"/></svg>
        LLMs & AI Validation Checkpoints
      </div>
      <p class="research-desc">
        Implementing rigorous checkpoints to audit AI-generated pipelines, mitigating silent bugs and batch artifacts, and assessing explainability maps (XAI) for biological reproducibility.
      </p>
    </div>
  </div>

  <!-- Software Section -->
  <div class="section-title">
    <span class="section-title-icon">
      <svg width="24" height="24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-code"><polyline points="16 18 22 12 16 6"></polyline><polyline points="8 6 2 12 8 18"></polyline></svg>
    </span>
    Developed Software & Systems
  </div>

  <div class="software-grid">
    
    <!-- ScanNeo -->
    <div class="software-card">
      <div class="card-top">
        <div class="card-header">
          <h3 class="card-title"><a href="https://github.com/ylab-hi/ScanNeo" target="_blank" rel="noopener noreferrer">ScanNeo</a></h3>
          <a href="https://github.com/ylab-hi/ScanNeo" class="card-github-icon" target="_blank" rel="noopener noreferrer" aria-label="ScanNeo GitHub">
            <svg width="20" height="20" fill="currentColor" viewBox="0 0 24 24"><path d="M12 2A10 10 0 0 0 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.9-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.9 1.52 2.34 1.07 2.91.83.1-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0 0 12 2z"/></svg>
          </a>
        </div>
        <p class="card-desc">Computational pipeline to identify insertions/deletions (indels) from RNA-Seq data and predict patient-specific, indel-derived neoantigens for personalized immunotherapy.</p>
        <div class="tag-container">
          <span class="tag primary">Bioinformatics</span>
          <span class="tag">RNA-Seq</span>
          <span class="tag">Neoantigen</span>
          <span class="tag">Python</span>
        </div>
      </div>
      <a href="https://academic.oup.com/bioinformatics/article/35/20/4159/5382215" class="card-action-link" target="_blank" rel="noopener noreferrer">
        Read Publication &rarr;
      </a>
    </div>

    <!-- ScanITD -->
    <div class="software-card">
      <div class="card-top">
        <div class="card-header">
          <h3 class="card-title"><a href="https://github.com/ylab-hi/ScanITD" target="_blank" rel="noopener noreferrer">ScanITD</a></h3>
          <a href="https://github.com/ylab-hi/ScanITD" class="card-github-icon" target="_blank" rel="noopener noreferrer" aria-label="ScanITD GitHub">
            <svg width="20" height="20" fill="currentColor" viewBox="0 0 24 24"><path d="M12 2A10 10 0 0 0 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.9-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.9 1.52 2.34 1.07 2.91.83.1-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0 0 12 2z"/></svg>
          </a>
        </div>
        <p class="card-desc">Algorithmic framework for detecting internal tandem duplications (ITD) from sequencing data with robust variant allele frequency estimation.</p>
        <div class="tag-container">
          <span class="tag primary">Genomics</span>
          <span class="tag">ITD Variant</span>
          <span class="tag">WGS/WES</span>
          <span class="tag">Python</span>
        </div>
      </div>
      <a href="https://academic.oup.com/gigascience/article/9/8/giaa089/5898622" class="card-action-link" target="_blank" rel="noopener noreferrer">
        Read Publication &rarr;
      </a>
    </div>

    <!-- ScanExitron -->
    <div class="software-card">
      <div class="card-top">
        <div class="card-header">
          <h3 class="card-title"><a href="https://github.com/ylab-hi/ScanExitron" target="_blank" rel="noopener noreferrer">ScanExitron</a></h3>
          <a href="https://github.com/ylab-hi/ScanExitron" class="card-github-icon" target="_blank" rel="noopener noreferrer" aria-label="ScanExitron GitHub">
            <svg width="20" height="20" fill="currentColor" viewBox="0 0 24 24"><path d="M12 2A10 10 0 0 0 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.9-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.9 1.52 2.34 1.07 2.91.83.1-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0 0 12 2z"/></svg>
          </a>
        </div>
        <p class="card-desc">An analytics toolkit for detecting and profiling exitron splicing (internally retained exonic region) events and predicting exitron-derived neoantigens.</p>
        <div class="tag-container">
          <span class="tag primary">Transcriptomics</span>
          <span class="tag">Exitrons</span>
          <span class="tag">Splicing</span>
          <span class="tag">R / Python</span>
        </div>
      </div>
      <a href="https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00297-7" class="card-action-link" target="_blank" rel="noopener noreferrer">
        Read Publication &rarr;
      </a>
    </div>

    <!-- OctopuSV & TentacleSV -->
    <div class="software-card">
      <div class="card-top">
        <div class="card-header">
          <h3 class="card-title"><a href="https://github.com/ylab-hi/OctopuSV" target="_blank" rel="noopener noreferrer">OctopuSV & TentacleSV</a></h3>
          <a href="https://github.com/ylab-hi/OctopuSV" class="card-github-icon" target="_blank" rel="noopener noreferrer" aria-label="OctopuSV GitHub">
            <svg width="20" height="20" fill="currentColor" viewBox="0 0 24 24"><path d="M12 2A10 10 0 0 0 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.9-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.9 1.52 2.34 1.07 2.91.83.1-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0 0 12 2z"/></svg>
          </a>
        </div>
        <p class="card-desc">Computational suites designed for multi-sample structural variant comparison, variant curation, and integration across varying sequencing technologies.</p>
        <div class="tag-container">
          <span class="tag primary">Genomics</span>
          <span class="tag">Structural Variant</span>
          <span class="tag">Long-read</span>
          <span class="tag">C++</span>
        </div>
      </div>
      <a href="https://github.com/ylab-hi/OctopuSV" class="card-action-link" target="_blank" rel="noopener noreferrer">
        View Repository &rarr;
      </a>
    </div>

    <!-- lncGSEA -->
    <div class="software-card">
      <div class="card-top">
        <div class="card-header">
          <h3 class="card-title"><a href="https://github.com/ylab-hi/lncGSEA" target="_blank" rel="noopener noreferrer">lncGSEA</a></h3>
          <a href="https://github.com/ylab-hi/lncGSEA" class="card-github-icon" target="_blank" rel="noopener noreferrer" aria-label="lncGSEA GitHub">
            <svg width="20" height="20" fill="currentColor" viewBox="0 0 24 24"><path d="M12 2A10 10 0 0 0 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.9-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.9 1.52 2.34 1.07 2.91.83.1-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0 0 12 2z"/></svg>
          </a>
        </div>
        <p class="card-desc">An R package to identify biological pathways and molecular mechanisms associated with long non-coding RNAs (lncRNAs) in cancer cohorts.</p>
        <div class="tag-container">
          <span class="tag primary">Functional Genomics</span>
          <span class="tag">lncRNA</span>
          <span class="tag">GSEA</span>
          <span class="tag">R</span>
        </div>
      </div>
      <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07873-y" class="card-action-link" target="_blank" rel="noopener noreferrer">
        Read Publication &rarr;
      </a>
    </div>

    <!-- HLA-IMPUTER -->
    <div class="software-card">
      <div class="card-top">
        <div class="card-header">
          <h3 class="card-title"><a href="https://github.com/dolittle007/HLA-IMPUTER" target="_blank" rel="noopener noreferrer">HLA-IMPUTER</a></h3>
          <a href="https://github.com/dolittle007/HLA-IMPUTER" class="card-github-icon" target="_blank" rel="noopener noreferrer" aria-label="HLA-IMPUTER GitHub">
            <svg width="20" height="20" fill="currentColor" viewBox="0 0 24 24"><path d="M12 2A10 10 0 0 0 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.9-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.9 1.52 2.34 1.07 2.91.83.1-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0 0 12 2z"/></svg>
          </a>
        </div>
        <p class="card-desc">An easy-to-use web application and association model for HLA allele imputation and analysis using custom population-specific reference panels.</p>
        <div class="tag-container">
          <span class="tag primary">Immunology</span>
          <span class="tag">HLA Imputation</span>
          <span class="tag">Shiny App</span>
          <span class="tag">R</span>
        </div>
      </div>
      <a href="https://academic.oup.com/bioinformatics/article/35/7/1244/5090124" class="card-action-link" target="_blank" rel="noopener noreferrer">
        Read Publication &rarr;
      </a>
    </div>

  </div>

  <!-- Publications Section -->
  <div class="section-title">
    <span class="section-title-icon">
      <svg width="24" height="24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-book-open"><path d="M2 3h6a4 4 0 0 1 4 4v14a3 3 0 0 0-3-3H2z"></path><path d="M22 3h-6a4 4 0 0 0-4 4v14a3 3 0 0 1 3-3h7z"></path></svg>
    </span>
    Selected Publications
  </div>

  <div class="pub-list">
    
    <!-- Pub 1 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.nature.com/articles/ng.3510" target="_blank" rel="noopener noreferrer">Genome-wide association meta-analysis in Chinese and European individuals identifies ten new loci associated with systemic lupus erythematosus</a>
      </div>
      <div class="pub-authors">
        DL Morris, Y Sheng, Y Zhang, <strong>YF Wang</strong>, Z Zhu, P Tombleson, L Chen, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Nature Genetics, 2016</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">367 Citations</span>
          <a href="https://www.nature.com/articles/ng.3510" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 2 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.nature.com/articles/s41467-021-21049-y" target="_blank" rel="noopener noreferrer">Identification of 38 novel loci for systemic lupus erythematosus and genetic heterogeneity between ancestral groups</a>
      </div>
      <div class="pub-authors">
        <strong>YF Wang</strong>, Y Zhang, Z Lin, H Zhang, <strong>TY Wang</strong>, Y Cao, DL Morris, Y Sheng, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Nature Communications, 2021</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">326 Citations</span>
          <a href="https://www.nature.com/articles/s41467-021-21049-y" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 3 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00297-7" target="_blank" rel="noopener noreferrer">A pan-cancer transcriptome analysis of exitron splicing identifies novel cancer driver genes and neoepitopes</a>
      </div>
      <div class="pub-authors">
        <strong>TY Wang</strong>, Q Liu, Y Ren, SK Alam, L Wang, Z Zhu, LH Hoeppner, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Molecular Cell, 2021</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">78 Citations</span>
          <a href="https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00297-7" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 4 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://academic.oup.com/bioinformatics/article/35/20/4159/5382215" target="_blank" rel="noopener noreferrer">ScanNeo: identifying indel-derived neoantigens using RNA-Seq data</a>
      </div>
      <div class="pub-authors">
        <strong>TY Wang</strong>, L Wang, SK Alam, LH Hoeppner, R Yang
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Bioinformatics, 2019</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">59 Citations</span>
          <a href="https://academic.oup.com/bioinformatics/article/35/20/4159/5382215" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 5 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.nature.com/articles/s41467-021-26614-7" target="_blank" rel="noopener noreferrer">Opposing transcriptional programs of KLF5 and AR emerge during therapy for advanced prostate cancer</a>
      </div>
      <div class="pub-authors">
        M Che, A Chaturvedi, SA Munro, SP Pitzen, A Ling, W Zhang, J Mentzer, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Nature Communications, 2021</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">45 Citations</span>
          <a href="https://www.nature.com/articles/s41467-021-26614-7" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 6 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://academic.oup.com/hmg/article/24/1/274/2900742" target="_blank" rel="noopener noreferrer">Meta-analysis of GWAS on two Chinese populations followed by replication identifies novel genetic variants on the X chromosome associated with systemic lupus erythematosus</a>
      </div>
      <div class="pub-authors">
        Y Zhang, J Zhang, J Yang, Y Wang, L Zhang, X Zuo, L Sun, HF Pan, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Human Molecular Genetics, 2015</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">45 Citations</span>
          <a href="https://academic.oup.com/hmg/article/24/1/274/2900742" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 7 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://ard.bmj.com/content/77/7/1078" target="_blank" rel="noopener noreferrer">Identification of ST3AGL4, MFHAS1, CSNK2A2 and CD226 as loci associated with systemic lupus erythematosus (SLE) and evaluation of SLE genetics in drug repositioning</a>
      </div>
      <div class="pub-authors">
        <strong>YF Wang</strong>, Y Zhang, Z Zhu, <strong>TY Wang</strong>, DL Morris, JJ Shen, H Zhang, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Annals of the Rheumatic Diseases, 2018</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">42 Citations</span>
          <a href="https://ard.bmj.com/content/77/7/1078" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 8 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://pubs.rsc.org/en/content/articlelanding/2011/mb/c1mb05018a" target="_blank" rel="noopener noreferrer">A predicted protein–protein interaction network of the filamentous fungus Neurospora crassa</a>
      </div>
      <div class="pub-authors">
        <strong>TY Wang</strong>, F He, QW Hu, Z Zhang
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Molecular BioSystems, 2011</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">27 Citations</span>
          <a href="https://pubs.rsc.org/en/content/articlelanding/2011/mb/c1mb05018a" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 10 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.science.org/doi/10.1126/sciadv.adl2804" target="_blank" rel="noopener noreferrer">EZH2 directly methylates PARP1 and regulates its activity in cancer</a>
      </div>
      <div class="pub-authors">
        Q Meng, J Shen, Y Ren, Q Liu, R Wang, Q Li, W Jiang, Q Wang, Y Zhang, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Science Advances, 2024</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">19 Citations</span>
          <a href="https://www.science.org/doi/10.1126/sciadv.adl2804" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 11 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://academic.oup.com/gigascience/article/9/8/giaa089/5898622" target="_blank" rel="noopener noreferrer">ScanITD: detecting internal tandem duplication with robust variant allele frequency estimation</a>
      </div>
      <div class="pub-authors">
        <strong>TY Wang</strong>, R Yang
      </div>
      <div class="pub-footer">
        <span class="pub-journal">GigaScience, 2020</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">19 Citations</span>
          <a href="https://academic.oup.com/gigascience/article/9/8/giaa089/5898622" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 15 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.cell.com/molecular-therapy-family/nucleic-acids/fulltext/S2162-2531(20)30005-7" target="_blank" rel="noopener noreferrer">Identification of regulatory modules that stratify lupus disease mechanism through integrating multi-omics data</a>
      </div>
      <div class="pub-authors">
        <strong>TY Wang</strong>, <strong>YF Wang</strong>, Y Zhang, JJ Shen, M Guo, J Yang, YL Lau, W Yang
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Molecular Therapy - Nucleic Acids, 2020</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">14 Citations</span>
          <a href="https://www.cell.com/molecular-therapy-family/nucleic-acids/fulltext/S2162-2531(20)30005-7" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 16 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.cell.com/star-protocols/fulltext/S2666-1667(21)00247-9" target="_blank" rel="noopener noreferrer">Integrated protocol for exitron and exitron-derived neoantigen identification using human RNA-seq data with ScanExitron and ScanNeo</a>
      </div>
      <div class="pub-authors">
        <strong>TY Wang</strong>, R Yang
      </div>
      <div class="pub-footer">
        <span class="pub-journal">STAR Protocols, 2021</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">11 Citations</span>
          <a href="https://www.cell.com/star-protocols/fulltext/S2666-1667(21)00247-9" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

    <!-- Pub 17 -->
    <div class="pub-item">
      <div class="pub-title">
        <a href="https://www.science.org/doi/10.1126/sciadv.eadk6989" target="_blank" rel="noopener noreferrer">Androgen receptor–regulated lncRNA PRCAT71 promotes AR signaling through the interaction with KHSRP in prostate cancer</a>
      </div>
      <div class="pub-authors">
        Y Yang, <strong>TY Wang</strong>, Q Li, J Lu, Y Ren, AB Weiner, J Fry, Q Liu, C Yum, ...
      </div>
      <div class="pub-footer">
        <span class="pub-journal">Science Advances, 2025</span>
        <div class="pub-stats">
          <span class="pub-citation-badge">8 Citations</span>
          <a href="https://www.science.org/doi/10.1126/sciadv.eadk6989" class="pub-link" target="_blank" rel="noopener noreferrer">Publisher Link</a>
        </div>
      </div>
    </div>

  </div>

</div>
