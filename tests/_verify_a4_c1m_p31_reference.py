"""Small independent controls for the P3.1 reference certification boundary."""
import numpy as np
from _a4_c1m_p31_reference import AffineReference,self_check
from _a4_c1m_p31_certify import gate,metric

def main():
    results=self_check()
    field={'h':np.array([0.,100.,200.]),'coefficients':np.array([[0.,0.],[0.,0.],[1e8,1e8],[1e11,1.1e11]])}
    for name,broken,direction in (
        ('nonaffine_rejected',dict(field,coefficients=np.array([[1.,0.],[0.,0.],[1e8,1e8],[1e11,1.1e11]])),[.5,.866,0.]),
        ('inexact_reduction_rejected',dict(field,coefficients=np.array([[0.,0.],[0.,0.],[1e8,1e8],[1e11,1.1e11+1.]])),[.5,.866,0.]),
        ('downward_out_of_scope_rejected',field,[-.5,.866,0.])):
        try:AffineReference(broken,[6471.,0.,0.],direction,80.61638494312905,10)
        except ValueError:results[name]=True
        else:raise AssertionError(name)
    exact={'xyz':np.array([[0.,0.,0.],[100.,0.,0.]]),'group':np.array([0.,100.])}
    common_bias=exact['xyz']+np.array([0.,.003,0.])
    a=metric(common_bias,exact['group'],exact);b=metric(common_bias.copy(),exact['group'],exact)
    assert np.array_equal(common_bias,common_bias.copy()) and not gate(a) and not gate(b)
    results['zero_cross_error_common_bias_rejected']=True
    assert not gate(dict(a,rms_m=0.,max_m=0.,landing_m=0.,group_rel=3e-6))
    results['group_only_bias_rejected']=True
    print(results);print('PASS '+str(len(results))+' checks')
if __name__=='__main__':main()
