from rf_pipelines import rf_pipelines_c

def clipper2d(Df, Dt, nt_chunk, sigma=3, niter=1, iter_sigma=3):
    return rf_pipelines_c.make_intensity_clipper2d(Df, Dt, nt_chunk, sigma, niter, iter_sigma)
