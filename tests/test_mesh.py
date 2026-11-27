import numpy as np
from pmwDPR.mesh import is_gmi, swath_mesh

def test_gmi_name():
    assert is_gmi('1C.GPM.GMI.X.HDF5')
    assert not is_gmi('1C.NPP.ATMS.X.HDF5')

def test_mesh_no_scatter():
    class Ax:
        def pcolormesh(self,x,y,z,**kw):
            assert np.isfinite(x).all() and np.isfinite(y).all()
            assert kw['shading']=='nearest'
            return z
    x,y=np.meshgrid(np.arange(5.),np.arange(5.)); x[0,0]=np.nan
    z=swath_mesh(Ax(),x,y,np.ones_like(x))
    assert z.mask[0,0] and not z.mask[-1,-1]
