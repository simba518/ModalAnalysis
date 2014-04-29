This application is used to compute modal analysis for Prof. Jin. 

The original email is:

Please help to compute/visualize the MA of the following data:
https://www.dropbox.com/sh/whu3r5ftcv8x4b4/DlqCUGDuOz

The format description can be found below:

".off" files contain hex-mesh data,  ".mass" and ".stiffness" files contain
the corresponding mass matrix and stiffness matrix data, respectively.

The format of .mass and .stiffness are the same. Since both of them are
sparse matrices, only elements with value not equal to 0 are recorded, as
described bellow:
"
#Row #Column
i a x  ------ For a line in the file, 'x' means the element value
corresponding to the ith row and ath column of the whole matrix.
"

Since I write the program all by myself, there are probabilities that the
matrices are calculated wrong. Any problems, please let me know.

Be careful, the mass matrix is not lumped.
