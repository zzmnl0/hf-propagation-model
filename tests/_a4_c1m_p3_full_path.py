"""Natural full-path metric contract, distinct from a fixed-length local scope."""
import numpy as np

def compare_case(d,jobs,out):
    """Versioned extended entry: preserve the frozen v2 local implementation."""
    from _a4_c1m_p3_compare import compare_case as local,read_trace
    if jobs[0]['scope']=='first_100_km':return local(d,jobs,out)
    jobs=sorted(jobs,key=lambda j:j['level'])
    if len(jobs)!=3:raise ValueError('Missing full-path integration level')
    traces={side:[read_trace(d,j,side) for j in jobs] for side in ('ah','ph')}
    row={'case_id':jobs[0]['job_id'].rsplit('_t',1)[0],'role':jobs[0]['role'],'scope':jobs[0]['scope'],'jobs':[j['job_id'] for j in jobs]}
    return analyze(row,traces)

def metrics(a,b,step=.25):
    from _a4_c1m_p3_compare import at,query
    end=min(float(a['s'][-1]),float(b['s'][-1]))
    if end<=0:raise ValueError('Empty full path')
    av=at(a,query(step,end));bv=at(b,query(step,end));distance=np.linalg.norm(av['xyz']-bv['xyz'],axis=1)*1000.
    return {'rms_m':float(np.sqrt(np.mean(distance**2))),'max_m':float(np.max(distance)),
            'landing_m':float(np.linalg.norm(a['xyz'][-1]-b['xyz'][-1])*1000.),
            'group_rel':float(abs(a['group'][-1]-b['group'][-1])/max(abs(a['group'][-1]),1e-15))}

def convergence(traces):
    from _a4_c1m_p3_compare import LIMITS
    d12=metrics(traces[0],traces[1]);d23=metrics(traces[1],traces[2])
    return {'delta12':d12,'delta23':d23,'monotonic':{k:d23[k]<=d12[k] for k in LIMITS},'last_change_small':{k:d23[k]<=.2*LIMITS[k] for k in LIMITS},'passed':all(d23[k]<=d12[k] and d23[k]<=.2*LIMITS[k] for k in LIMITS),'roundoff_waiver_applied':False}

def analyze(row,traces):
    from _a4_c1m_p3_compare import LIMITS
    actual=metrics(traces['ah'][-1],traces['ph'][-1]);integration={s:convergence(ts) for s,ts in traces.items()}
    eligible=all(t['termination']=='ground' for t in traces['ah']) and all(t['termination']==1 for t in traces['ph'])
    row.update(actual=actual,integration=integration,topology={'expected':'ground','eligible':eligible,'matched':eligible,'coverage_complete':eligible},
               natural_endpoints={s:[{'xyz_km':t['xyz'][-1].tolist(),'s_km':float(t['s'][-1]),'group_km':float(t['group'][-1]),'termination':t['termination']} for t in ts] for s,ts in traces.items()},
               common_arc_km=min(traces['ah'][-1]['s'][-1],traces['ph'][-1]['s'][-1]),
               reconstruction={'validated':False,'reason':'Full-path independent reconstruction bound not established for this input; local P2 bounds are inapplicable'},
               checks={'thresholds':all(actual[k]<=LIMITS[k] for k in LIMITS),'integration':all(v['passed'] for v in integration.values()),'scope':eligible,'reconstruction':False},strict_passed=False)
    return row

def self_check():
    from _a4_c1m_p3_compare import query
    def curve(end):
        s=query(1.,end);xyz=np.c_[200*np.sin(s/200),200*(1-np.cos(s/200)),s*0];tangent=np.c_[np.cos(s/200),np.sin(s/200),s*0]
        return {'s':s,'xyz':xyz,'tangent':tangent,'p':tangent,'group':s,'termination':'ground'}
    a,b=curve(100.),curve(101.);m=metrics(a,b)
    assert m['rms_m']<1e-7 and m['landing_m']>999 and m['group_rel']==.01
    short=curve(90.);short['termination']='max_path';ph=curve(100.);ph['termination']=1
    rejected=analyze({}, {'ah':[a,a,short],'ph':[ph,ph,ph]})
    assert not rejected['topology']['eligible'] and not rejected['strict_passed']
    linear={**a,'xyz':np.c_[a['s'],a['s']*0,a['s']*0],'tangent':np.tile([1.,0,0],(len(a['s']),1)),'p':np.tile([1.,0,0],(len(a['s']),1))}
    assert metrics(linear,linear)==dict.fromkeys(('rms_m','max_m','landing_m','group_rel'),0.)
    return {'straight_exact':True,'curved_common_arc_separated_from_natural_endpoint':True,'group_uses_each_natural_end':True,'early_termination_rejected':True,'full_path_reference_validated':False}
