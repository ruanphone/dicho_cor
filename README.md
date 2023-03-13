# dicho_cor
## User’s manual for dicho_cor_x and dicho_cor_xy

Assume a continuous variable x to be a Gaussian mixture, R codes dicho_cor_x() calculates the optimal cut point h of x that maximize the point‐biserial correlation between x and the dichotomized x (denoted by xd(h)).  
To run the program, you need to install the following R packages first: mclust, plyr, huxtable.  
And then input the target X into the parenthesis (e.g., dicho_cor_x(X)).  
The cut point h, Max(Corr(xd (h),x) and the results from Gaussian mixture classifier will be outputted.  

Be noted that no missing value is allowed.

Further assume that outcome y is continuous and has piece‐wise linear relation with exposure x, the R codes dicho_cor_xy() calculates the optimal cut point h of x that maximize the point‐biserial correlation between y and the dichotomized x (xd(h)).  
To run the program, you need to install the following R packages first: mclust, plyr, huxtable.And the input X and Y (e.g., , dicho_cor_xy(X,Y)).  
The cut point h, Max(Corr(x_d (h),y)) and and the results from Gaussian mixture classifier will be outputted. No missing value is allowed.  

<https://www.dropbox.com/s/app3a06dpwdkuct/dicho_cor_x.txt>  
<https://www.dropbox.com/s/yn6137g3jbbiiy1/dicho_cor_xy.txt>  
The flowcharts of the algorithms are shown as followed.  

