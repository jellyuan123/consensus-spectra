# Consensus Spectra Generation and Evaluation

**Author:** Jasmine Schaber  
**Course:** CS166 – Intro to Computational Biology  
**Date:** December 2024  

## Overview  
This project investigates whether generating **consensus spectra** from multiple mass spectrometry samples can improve molecular **candidate ranking** accuracy using tools like MSGym, SKLearn, and JESTR.

## Goals  
- Implement consensus spectra generation using binning and clustering  
- Evaluate impact on candidate rankings via JESTR  
- Explore how factors like collision energy and sample size affect consensus quality  

## Methodology  

### Data Preprocessing  
- Used **MSGym** to download and prepare large spectral datasets.

### Consensus Generation  
- Binned m/z peaks using 0.02 m/z width.  
- Retained bins with peaks from >25% of the sample spectra.  
- Averaged intensities of retained peaks to form the consensus.  
- Extracted precursor m/z values; averaged or selected the majority value based on variance.

### Clustering  
- Employed **SKLearn** to cluster peaks by setting the cluster count to match the average number of peaks per spectrum.

### Candidate Ranking  
- Selected 10 random consensus spectra.  
- Converted to `.mgf` format.  
- Submitted to **JESTR** for SMILES candidate ranking evaluation.

## Results  
- 7 consensus spectra sets successfully ranked.  
- Only one spectrum achieved a top-10 rank (rank = 3).  
- Others ranked poorly (ranked 25 or in the hundreds).  
- ~32,000 spectra generated from ~36,000; causes of failures unknown but not data-limiting.

## Challenges  
- Generally poor ranking performance.  
- Unexpected cluster ranges and potentially flawed binning.  
- Limited sample tracking (size per molecule not recorded).  
- JESTR access was limited to peer support.  
- Full generation took over a day on standard hardware.

## Future Work  
- Retrain **JESTR** using consensus spectra to assess learning improvements.  
- Refine binning and clustering logic.  
- Track sample sizes during spectra generation.  
- Expand the number of test trials and automate logging/error capture.

## References  
- Roman Bushuiev et al., *MassSpecGym: A benchmark for the discovery and identification of molecules*, 2024.  
- Xiyang Luo et al., *A Comprehensive Evaluation of Consensus Spectrum Generation Methods in Proteomics*, J. Proteome Res., 2022, 21(6), 1566–1574. DOI: [10.1021/acs.jproteome.2c
