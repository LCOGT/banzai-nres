import numpy as np

import logging

from banzai.stages import Stage
from banzai.utils import qc

logger = logging.getLogger('banzai')

class CheckIfFluxInTraces(Stage):
    """
    Verify that most of the flux falls within the traces

    @author:mjohnson
    """

    def __init__(self,runtime_context):
        super(CheckIfFluxInTraces, self).__init__(runtime_context)
    
    def do_stage(self,image):
        self.check_flux_in_and_out_of_traces(image, flux_ratio)
        self.check_flux_by_trace(image, flux_by_trace)
        qc_results = {'flux_ratio':flux_ratio, 'flux_by_trace':flux_by_trace}
        qc.save_qc_results(self.runtime_context, qc_results, image)
        return image

    def check_flux_in_and_out_of_traces(self,image):
        traces = np.copy(image.traces)
        traces[traces != 0] = 1
        flux_in_traces = np.sum(image.data*traces)
        flux_outside_traces = np.sum(image.data*(traces-1)*(-1.))
        flux_ratio = flux_in_traces/flux_outside_traces
        return flux_ratio

    def check_flux_by_trace(self,image):
        trace_ids = np.arange(1, image.num_traces + 1)
        flux_by_trace = np.ones_like(trace_ids)
        for i, trace_id in enumerate(trace_ids): 
            flux_by_trace[i] = np.sum(image.data[image.traces == i])
        return flux_by_trace
