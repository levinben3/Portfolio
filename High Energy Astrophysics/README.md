## Astrophysics Projects
All projects related to work with the Columbia High Energy Astrophysics Laboratory is housed in this repository. The Jupyter notebooks python scripts in this directory are all related to a key project studying the ESA XMM-Newton Satellite. The project aims to stack energy spectra of thousands of faint sources in the galactic center all housed in the XMM-Newton Serendipitous Source Catalog, which compiles all detections from 2000-2023. The stacking requires a variety of scripts, most notably including a novel grouping pipeline that will be key in final analysis. The project aims to compelete a never-accomplished characteriztion of all relevant sources in this archive. 

Project and Paper Links Below

| File            | Explanation                                                                |
| ----------------- | ------------------------------------------------------------------ |
| download_script | Web script to download, label, and sort spectral files from Nasa's HEASARC XAMIN archive. |
| merge_pipeline | Novel grouping pipeline employing contour map grouping, sigma clipping and k-means-clustering-like scripts, and other grouping based on source properties |
| plots | A script that works in tandem with the main pipeline, this provides many interesting and robust visualizations of grouping results to both test and document the pipeline. |
| quantile_analysis | A test script on a novel replacement for the commonly-used hardness ratio, which describes how energetic a source is. This is a much more robust way to characterize source, while allowing for error calculation and fewer individual photon detections. |

##Project Links
* Project Presentation and Guide: [![project_presentation]](https://docs.google.com/presentation/d/1HL4Q8onpGobTwm2RAeGiQ0oTclORFduwX6KvxPzGspk/edit?usp=sharing)
* Astronomer's Telegram Discovery Announcement: [![atel](https://www.astronomerstelegram.org/?read=17087)
