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


def test_text_catalog():
    from pmwmatch.automatic import discover, interval, orbit
    class Response:
        status_code=200
        text='/gpmdata/2015/01/01/radar/2A.GPM.DPR.V10-20260310.20150101-S012932-E030203.004783.V08A.nc\n'
        def raise_for_status(self): pass
    class Session:
        def get(self,url,timeout):
            assert '/text/gpmdata/' in url
            return Response()
    a,b=interval(Response.text)
    urls=discover(Session(),a,b,'radar',lambda n:True)
    assert len(urls)==1 and orbit(urls[0])==4783 and '/text/' not in urls[0]


def test_track_events_dateline():
    import numpy as np
    from pmwmatch.trackplots import events
    a=(np.array([0.]),np.array([179.99]),np.array([100.]),'track')
    b=(np.array([0.]),np.array([-179.99]),np.array([110.]),'track')
    assert events(a,b,7.5,15)==[(0,0)]
    assert events(a,b,0.1,15)==[]
