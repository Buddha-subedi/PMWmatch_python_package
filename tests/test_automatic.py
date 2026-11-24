import numpy as np
from pmwmatch.automatic import map_group,interval

def test_temporal_candidate_and_nan():
    s=dict(lat=np.array([[0.,0.]]),lon=np.array([[0.,.01]]),time=np.array([[1000.,60.]]),data=np.array([[[200.,210.],[250.,np.nan]]]))
    t=dict(lat=np.array([[0.,10.]]),lon=np.array([[0.,0.]]),time=np.array([[0.,0.]]))
    r=map_group(s,t,1,15)
    assert r['pmw_tb'][0,0,0]==250
    assert np.isnan(r['pmw_tb'][0,0,1])
    assert np.isnan(r['pmw_tb'][0,1]).all()
    assert r['delta_t_minutes'][0,0,0]==1

def test_midnight():
    a,b=interval('1C.F16.SSMIS.XCAL2021.20230101-S235900-E000100.123.V07A.HDF5')
    assert (b-a).total_seconds()==120
