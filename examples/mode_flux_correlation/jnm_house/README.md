## Computing time integral of correlation.

First compile the codes:

    g++ cor_fft.cpp -lfftw3 -o calc
    g++ integrate.cpp -o integrate
    g++ integrate_ave.cpp -o integrate_ave #Same as the above code, except has an option for running average to smooth the time integrating curve
    
Use FFT to calculate the correlation as a function of time:

    ./calc
    
Which makes cordata.txt. We now time-integrate this with

    ./integrate ### JUST DO ./integrate_ave SINCE IT MAKES BOTH NORMAL AND RUNNING AVE PLOTS.

which makes the time-integrated GK expression in area_evol_CCs.txt. A running average (smoother) version of this is done by:

    ./integrate_ave
    
To plot the time-integrated correlation function, do:

    gnuplot < gnuplot_script
    
To plot the running average version, do:

    gnuplot < gnuplot_script_ave

All together:

    ./calc
    ./integrate_ave
    gnuplot < gnuplot_script
    gnuplot < gnuplot_script_ave
