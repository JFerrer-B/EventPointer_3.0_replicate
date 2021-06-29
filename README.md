# EventPointer 3.0

[![platforms](http://bioconductor.org/shields/availability/3.8/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)
[![build](http://bioconductor.org/shields/build/release/bioc/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)
[![updated](http://bioconductor.org/shields/lastcommit/release/bioc/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)
[![yearsBioC](http://bioconductor.org/shields/years-in-bioc/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)


Code to replicate results depicted in the article:

**EventPointer3.0: flexible and accurate splicing analysisthat includes studying the differential usage of protein-domains**

*EventPointer* is an R package to identify alternative splicing events 
		that involve either simple (case-control experiment) or complex experimental designs 
		such as time course experiments and studies including paired-samples. The algorithm can
		be used to analyze data from either junction arrays (Affymetrix Arrays) or sequencing data (RNA-Seq). 
		The software returns a data.frame with the detected alternative splicing 
		events: gene name, type of event (cassette, alternative 3',...,etc), genomic 
		position, statistical significance and increment of the percent spliced in (Delta PSI) for all 
		the events.
		The algorithm can generate a series of files to visualize the detected alternative 
		splicing events in IGV. This eases the interpretation of results and the design 
		of primers for standard PCR validation.
		
This repository has three main directories:

* simulation_data: code to replicate the analysis performed with simulated data.
* HVS_dataset: code to replicate the analysis performed with the HVS dataset.
* cx4945_dataset: code to replicate the analysis performed with the cx4945dataset

