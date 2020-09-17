### Objective: To create a fully autonomous data analysis pipeline to look for Fast Radio Bursts (FRBs) working 24*7 for all observing frequency bands of Giant Metrerwave Radio Telescope. 

### Pipeline flow : 171 sec data split, SHM Blocks and design, different frequency bands/time resolution and conversion to 1k@1.3ms. 

SHM Design: The telescope has a buffer of 2GB. We have created a SHM of 

Overlap: How it is calculated? 

Python file: What does it do? Matched filtering

Overallflow of pipeline: Reading of telescope metadata, Reading from SHM in 8 blocks(nsamples??), sending to GPU, output peaks files,
                         Python match filtering, CSV file output, FETCH pipeline(CNNs)
