import argparse,json
from pathlib import Path
import yaml
from .readers import inspect_swaths,read_pmw,read_dpr
from .workflows import PMWDPR
from .mesh import is_gmi
from .visualize import plot_collocation,plot_orbit_tracks
from .orbits import find_intersections
from .download import download_urls,download_pps

def main():
    p=argparse.ArgumentParser(prog='pmwDPR')
    sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser("auto",help="Interactive PPS discovery and matching")
    a=sub.add_parser('inspect');a.add_argument('file')
    a=sub.add_parser('run');a.add_argument('config')
    a=sub.add_parser('plot');a.add_argument('npz');a.add_argument('--channel',type=int,default=0);a.add_argument('--step',type=int,default=1)
    a=sub.add_parser('download');a.add_argument('urls_file');a.add_argument('--out',required=True);a.add_argument('--pps',action='store_true')
    args=p.parse_args()
    if args.cmd=='auto':
        print(PMWDPR())
    elif args.cmd=='inspect':
        for name,shape,dtype in inspect_swaths(args.file):print(f'{name:65} {str(shape):22} {dtype}')
    elif args.cmd=='plot':print(plot_collocation(args.npz,channel=args.channel,step=args.step))
    elif args.cmd=='download':
        urls=[x.strip() for x in Path(args.urls_file).read_text().splitlines() if x.strip() and not x.startswith('#')]
        print('\n'.join((download_pps if args.pps else download_urls)(urls,args.out)))
    else:
        cfg=yaml.safe_load(Path(args.config).read_text())
        kind=cfg['workflow'].upper();kw=cfg.get('matching',{})
        gmi=is_gmi(cfg['pmw_file']) or str(cfg.get('sensor','')).upper()=='GMI'
        if gmi: kw['same_platform']=True
        if kind=='PMWDPR':
            result=PMWDPR(cfg['pmw_file'],cfg['reference_file'],cfg['pmw_reader'],cfg['reference_reader'],cfg['output_file'],**kw)
            if not result.get('skipped') and not gmi and cfg.get('orbit_figure'):
                pmw=read_pmw(cfg['pmw_file'],cfg['pmw_reader']);dpr=read_dpr(cfg['reference_file'],cfg['reference_reader'])
                crossings=[] if kw.get('same_platform',False) else find_intersections(pmw,dpr,kw.get('coincidence_minutes',7))
                plot_orbit_tracks(pmw,dpr,crossings,cfg['orbit_figure'])
        else:raise ValueError('workflow must be PMWDPR')
        print(json.dumps(result,indent=2))
        if not result.get('skipped') and not gmi and cfg.get('plot',True):print('Figure:',plot_collocation(cfg['output_file']))
if __name__=='__main__':main()
