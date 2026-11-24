import numpy as np
from pmwDPR.matching import map_nearest

def test_no_temporal_interpolation():
    src=dict(lat=np.array([[0.,0.]]),lon=np.array([[0.,1.]]),time=np.array([[0.,120.]]),data=np.array([[[200.],[250.]]]))
    tgt=dict(lat=np.array([[0.,0.]]),lon=np.array([[0.,1.]]),time=np.array([[0.,0.]]))
    r=map_nearest(src,tgt,max_distance_km=214,pixel_time_minutes=1)
    assert r['tb'][0,0,0]==200
    assert np.isnan(r['tb'][0,1,0])

def test_distance_cutoff():
    src=dict(lat=np.array([[0.]]),lon=np.array([[0.]]),time=np.array([[0.]]),data=np.array([[[200.]]]))
    tgt=dict(lat=np.array([[10.]]),lon=np.array([[0.]]),time=np.array([[0.]]))
    r=map_nearest(src,tgt)
    assert np.isnan(r['tb'][0,0,0])
