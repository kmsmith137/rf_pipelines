Some day, there will be systematic documentation here!

In the meantime, here are some notes on transforms which are new,
or have been changed recently:

 - `spline_detrender(nt_chunk, axis, nbins, epsilon = 3.0e-4)`

   This is an experimental transform designed to address the high
   computational cost of the polynomial_detrender.
   
   A spline_detrender with N bins should be roughly equivalent to a
   polynomial_detrender with degree (2N+2), but its computational cost
   should be independent of N.
    
   The 'epsilon' parameter regulates the spline fit by penalizing large time derivatives.
   If epsilon is too small, then overfitting may occur in regions with sparse weights.
   If epsilon is too large, then the fitter may have difficulty "keeping up" with rapid
   variations in the data.  I think that 3.0e-4 (the default) is a reasonable choice of
   epsilon, but I haven't experimented systematically.
    
   FIXME: currently, the only allowed axis type is AXIS_FREQ (=1).
   I propose that we experiment with replacing the polynomial_detrender(AXIS_FREQ),
   and if this looks good, then I'll implement AXIS_TIME and AXIS_NONE.
