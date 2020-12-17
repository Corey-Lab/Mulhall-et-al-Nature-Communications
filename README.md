# Mulhall-et-al-Nature-Communications

This is the analysis/simulation code used in the paper:

Single-molecule Force Spectroscopy Reveals the Dynamic Strength of the Hair-Cell Tip-Link Connection 

All the scripts have been exectued in Matlab 2017 and 2018. Fitting functions require the curve fitting toolbox and the optimization toolbox. To create probabiity density function for rupture forces using the ksdensity function requires the statistics and machine learning toolbox.

The folder DataAnalysisExport has all the functions that can be used to analyze rupture force spectroscopy for a variety of different models. The file DataAnalysisScript.m has various sections that can be used to analyze the data from our paper. There also sections in the script where you can import the rupture force data from the available excel spreadsheet. There are sections to import the monomer and dimer data. Once the raw rupture force data is imported the script as sections to calculate the mostprobable rupture force. This function calculates most probable rupture a few different ways: making a standard histogram, which uses Freedman-Diaconis rule for binwidth, using the cumulative rupture force distribution to estimate most probable rupture force, using the ksdensity function in matlab, using the average shifted histogram method, and also calculating the median rupture force. The error in force outputed is either the histogram binwidth, the kernel smoothing width parameter in ksdensity, or the median aboslute deviation. After this you define your initial guess parameters for fitting, then fit the data. NOTE: the script is not really ment to be run as a hole, but rather each section exectued individually. Depending on what data you are trying fit and what model you want to use. You can copy and past the monomer and dimer data from figure 3 in the paper and run the appropriate sections to obtain the plots in the paper.


The SimulationExport folder has all the functions needed to run Monte-Carlo simulations of the oscilitory stimulation. The SimulationScript.m file is a script that will run a batch of simulations. The scriptfile is fairly simple and should be able to be ran by looking at the comments. The output of the simulation will be 3D matrix of lifetime values for each simulation. Each column is the lifetime for each replicate for a given amplitude. Each z index is for a given frequency.

All code was written by Andrew Ward. If you have any questions on using the code please email andyrward@gmail.com.
