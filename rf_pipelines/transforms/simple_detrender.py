from rf_pipelines import rf_pipelines_c

def simple_detrender(nt_detrend, axis=1, deg=0, epsilon=1.0e-2):
    if axis == 1:
        return rf_pipelines_c.make_polynomial_detrender_time_axis(nt_detrend, deg, epsilon)
    elif axis == 0:
        return rf_pipelines_c.make_polynomial_detrender_freq_axis(nt_detrend, deg, epsilon)
    else:
        raise RuntimeError("axis must be 0 (along freq) or 1 (along time).")
