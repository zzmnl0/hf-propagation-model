"""P3 pointwise state diagnostics and delivery plots; never changes raw results."""
from pathlib import Path
import argparse,csv,json,shutil
import numpy as np
from scipy.interpolate import CubicHermiteSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import _a4_c1m_p3_compare as c

def artifacts(run,analysis='p3_v2',version='p3_artifacts_v1'):
    d=run/'experiments';out=run/'analysis'/version
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);shutil.copy2(Path(__file__),out/'execution_source.py');rows=[]
    for j in c.exp.read(d/'manifest.json')['jobs']:
        g=c.common.mat(Path(j['grid_mat']));h=np.asarray(g['height_km']);ne=np.asarray(g['ne_profile_m3']);spline=CubicHermiteSpline(h,ne,np.gradient(ne,h,edge_order=2))
        for side in ('ah','ph'):
            if side=='ph' and j['role']=='R5_fd':continue
            try:t=c.read_trace(d,j,side)
            except (ValueError,FileNotFoundError):continue
            heights=np.linalg.norm(t['xyz'],axis=1)-6371.;input_ne=spline(np.clip(heights,h[0],h[-1]));input_ne[(heights<h[0])|(heights>h[-1])]=0.
            p2=np.sum(t['p']**2,axis=1);residual=p2-(1-c.exp.COEFFICIENT*input_ne/(j['frequency_MHz']*1e6)**2)
            angular=np.cross(t['xyz'],t['p']);angular_drift=np.linalg.norm(angular-angular[0],axis=1)/max(np.linalg.norm(angular[0]),1e-30)
            output_ne=np.full(len(heights),np.nan);n2=np.full(len(heights),np.nan)
            if side=='ph':
                raw=c.common.mat(Path(t['raw_path']));path=raw['ray_path_data'];output_ne=np.asarray(path['electron_density']).reshape(-1)*1e6;n2=np.asarray(path['refractive_index']).reshape(-1)**2
            else:output_ne=input_ne;n2=1-c.exp.COEFFICIENT*input_ne/(j['frequency_MHz']*1e6)**2
            pp=out/(j['job_id']+'_'+side+'_native.csv')
            with pp.open('w',newline='',encoding='utf-8') as f:
                w=csv.writer(f);w.writerow(['s_km','x_km','y_km','z_km','p_x','p_y','p_z','height_km','group_km','Ne_input_signed_m3','Ne_output_m3','p2','n2_output','dispersion_residual_input','angular_momentum_relative_drift'])
                for i in range(len(heights)):w.writerow([t['s'][i],*t['xyz'][i],*t['p'][i],heights[i],t['group'][i],input_ne[i],output_ne[i],p2[i],n2[i],residual[i],angular_drift[i]])
            use=t['s']<=100.+1e-8 if j['scope']=='first_100_km' else np.ones(len(heights),dtype=bool)
            rows.append({'job_id':j['job_id'],'side':side,'scope':j['scope'],'native_points':len(heights),'residual_max':float(np.max(np.abs(residual[use]))),'angular_relative_drift_max':float(np.max(angular_drift[use])),'Ne_difference_max_m3':float(np.max(np.abs(output_ne[use]-input_ne[use]))),'raw_sha256':t['raw_sha256']})
    c.p2.write_csv(out/'native_summary.csv',rows)
    summary=run/'analysis'/analysis/'summary.json'
    if summary.exists():
        cases=c.exp.read(summary)['cases'];computed=[r for r in cases if 'actual' in r]
        c.p2.write_csv(out/'case_metrics.csv',[{'case_id':r['case_id'],'role':r['role'],'scope':r['scope'],'passed':r['strict_passed'],**r.get('actual',{}),'error':r.get('error','')} for r in cases])
        fig,ax=plt.subplots(figsize=(12,5));x=np.arange(len(computed));ax.semilogy(x,[r['actual']['rms_m'] for r in computed],'o-',label='RMS');ax.semilogy(x,[r['actual']['max_m'] for r in computed],'s-',label='Maximum');ax.axhline(10,color='k',ls='--',label='RMS limit');ax.set_xticks(x,[r['case_id'] for r in computed],rotation=90,fontsize=6);ax.set_ylabel('Position discrepancy (m)');ax.legend();fig.tight_layout();fig.savefig(out/'position_discrepancies.png',dpi=150);plt.close(fig)
    scans=c.exp.read(run/'analysis/p3_r5_v1/summary.json');fig,axs=plt.subplots(1,2,figsize=(10,4))
    for r in scans['fd']:
        if 'controls' not in r:continue
        axs[0].plot([.25,.125],[r['controls']['delta12']['rms_m'],r['controls']['delta23']['rms_m']],'o-',label=r['base_job_id'])
    axs[0].set(xlabel='Finer FD step (km)',ylabel='Successive RMS change (m)',title='FD controls');axs[0].legend(fontsize=6)
    for r in scans['grid']:
        if 'successive_changes' not in r:continue
        axs[1].plot([1,.5,.25],[x['rms_m'] for x in r['successive_changes']],'o-',label=str(r['case']))
    axs[1].set(xlabel='Finer height grid (km)',ylabel='Successive RMS change (m)',title='Affine grid controls');axs[1].legend(fontsize=5);fig.tight_layout();fig.savefig(out/'single_factor_controls.png',dpi=150);plt.close(fig)
    c.exp.write(out/'manifest.json',{'source_sha256':c.common.sha(Path(__file__)),'files':{p.name:c.common.sha(p) for p in out.iterdir() if p.is_file()},'raw_modified':False})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--analysis',default='p3_v2');p.add_argument('--version',default='p3_artifacts_v1');a=p.parse_args();artifacts(a.run.resolve(),a.analysis,a.version)
