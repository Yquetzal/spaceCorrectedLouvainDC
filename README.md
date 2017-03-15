# spaceCorrectedLouvainDC
Doubly constrained version of the space corrected community detection.

Check the ExampleUsageSpatialNullM to understand how to use it.

When calling the function that compute the null model, it can be useful to change "plot" and "printDebug" to better understand the process.

An important parameter is "roundDecimal". It corresponds to the size of the bins when computing the deterrence function. Necessary to change it to adapt to the distance format !
If you want to check if the deterrence function looks good or not, you can set the "plot" argument at true, this will display the deterrence function during first and last steps.
