import sys
import rf_pipelines
import FDMT
import matplotlib.pyplot as plt
import cart
class fdmt_transform(rf_pipelines.py_wi_transform):

    def __init__(self):

        rf_pipelines.py_wi_transform.__init__(self)


    def set_stream(self, stream):
        pass


    def start_substream(self, isubstream, t0):
        pass

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        I = intensity.copy()
        D = FDMT.FDMT(I, 400, 800, 64*1024, 'int64')

        print cart.argmaxnd(D)

    def end_substream(self):
        pass
