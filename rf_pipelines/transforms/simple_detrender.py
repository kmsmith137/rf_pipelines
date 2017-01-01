# A "simple detrender" is a time-axis polynomial fitter with degree zero.
# This file will be removed soon, in favor of calling make_polynomial_detrender_time_axis() directly.

from rf_pipelines import rf_pipelines_c

def simple_detrender(nt_detrend):
    print >>sys.stderr, 'make_simple_detrender(): this function is now deprecated'
    return rf_pipelines_c.make_polynomial_detrender_time_axis(nt_detrend, 0, 1.0e-2)
