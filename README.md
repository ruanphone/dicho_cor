# dicho_cor
## **User’s manual for dicho_cor_x and dicho_cor_xy**

Assume a continuous variable x to be a Gaussian mixture, R codes _**dicho_cor_x()**_ calculates the optimal cut point h of x that maximize the point‐biserial correlation between x and the dichotomized x (denoted by xd(h)).  
To run the program, you need to install the following R packages first: mclust, plyr, huxtable.  
And then input the target X into the parenthesis (e.g., _**dicho_cor_x(X)**_).  
The cut point h, Max(Corr(xd (h),x) and the results from Gaussian mixture classifier will be outputted.  

Be noted that no missing value is allowed.

Further assume that outcome y is continuous and has piece‐wise linear relation with exposure x, the R codes _**dicho_cor_xy()**_ calculates the optimal cut point h of x that maximize the point‐biserial correlation between y and the dichotomized x (xd(h)).  
To run the program, you need to install the following R packages first: mclust, plyr, huxtable.And the input X and Y (e.g., , _**dicho_cor_xy(X,Y)**_).  
The cut point h, Max(Corr(x_d (h),y)) and and the results from Gaussian mixture classifier will be outputted. No missing value is allowed.  

<https://www.dropbox.com/s/app3a06dpwdkuct/dicho_cor_x.txt>  
<https://www.dropbox.com/s/yn6137g3jbbiiy1/dicho_cor_xy.txt>  
The flowcharts of the algorithms are shown as followed.  

The three data sets used in the paper is available through:  
[Data 1 HbA1c vs glucose](https://www.dropbox.com/s/rz4hg2q9zy3kk40/data1.csv?dl=0)  
[Data 2. 921 earthquake](https://www.dropbox.com/s/z7dud1o5f35zl9z/921report2.htm)  
[Data 3. Height and Weight](https://www.dropbox.com/s/kipgp1hmigorqdx/weight_height.csv)  
[Data 4 Insurance charge](https://data.world/bob-wakefield/insurance)  

![image](https://github.com/iblian/dicho_cor/blob/main/Image.png)
