---
title: 'fastq-dl: efficiently download sequences from ENA and SRA'
tags:
  - fastq
  - download
  - python
  - bioinformatics
authors:
  - name: Robert A. Petit III
    orcid: 0000-0002-1350-9426
    affiliation: "1, 2"
  - name: Michael B. Hall
    orcid: 0000-0003-3683-6208
    affiliation: "3"
  - name: Gerry Tonkin-Hill
    orcid: 0000-0002-1350-9426
    affiliation: "4"
  - name: Jie Zhu
    affiliation: "5"
  - name: Timothy D. Read
    orcid: 0000-0001-8966-9680
    affiliation: "2"
affiliations:
  - name: Wyoming Public Health Laboratory, Wyoming Department of Health, Cheyenne, Wyoming, USA
    index: 1
  - name: Division of Infectious Diseases, Department of Medicine, Emory University School of Medicine, Atlanta, Georgia, USA
    index: 2
  - name: Department of Microbiology and Immunology, Peter Doherty Institute for Infection and Immunity, The University of Melbourne, Melbourne, Australia
    index: 3
  - name: Department of Biostatistics, University of Oslo, Oslo, Norway
    index: 4
  - name: Li Ka Shing Institute of Health Sciences, Faculty of Medicine, The Chinese University of Hong Kong, Hong Kong SAR, PR China
    index: 5
date: 17 June 2023
bibliography: paper.bib
---

# Summary

High-throughput sequencing technologies have revolutionized the field of genomics, enabling
researchers to generate vast amounts of data quickly and at relatively low cost. The European
Nucleotide Archive (ENA) [@Burgin_2023] and the Sequence Read Archive (SRA) [@Katz_2022] are
two major repositories for publicly hosting next-generation sequencing data from many research
projects. Retrieving sequences from these repositories is often a multi-step process and
difficult for researchers who lack experience with bioinformatics. fastq-dl is a bioinformatic
tool that simplifies the process of downloading sequences from SRA and ENA.

fastq-dl is written in Python and is designed to be user-friendly and simple to use. Users can
submit queries to the ENA, via a REST API [@Burgin_2023], or SRA, via pysradb [@Choudhary_2019],
with fallback mechanisms in the event either repository is down. fastq-dl supports a range of
query types, including taxon ids, species names, and accessions, including BioSample, BioProject,
Experiment, and Run Accessions. A query will return metadata for each hit and save this metadata
to a tab-delimited file. Unless disabled by the user, fastq-dl will then proceed to download
available sequences for each hit of the query. If using ENA, raw FASTQs are downloaded using
their available FTP service, otherwise fasterq-dump [@Katz_2022] is used to download from SRA.
In the event a repository is unresponsive, download attempts will be made against the other
repository. When an Experiment or BioSample has multiple Run accessions associated with it,
users can optionally choose to merge these Run accessions. Upon completion, users are provided
with a summary file, a metadata file and FASTQ files per-query hit.

fastq-dl is a convenient bioinformatic tool that simplifies the process of retrieving FASTQ files
from ENA and SRA. It was developed to be easy to use and accessible to researchers from all
backgrounds. By facilitating efficient downloading of publicly available FASTQ files, users can
easily integrate these data into their own research. fastq-dl is available from PyPI and Bioconda
for simple installation, and the source code is available at https://github.com/rpetit3/fastq-dl.  

# Funding

This project was partially supported by the Georgia Emerging Infections Program and the Wyoming Department of Health

# References
