"""Explicit HTTPS downloader. Source catalog URLs must be verified by the user."""
from pathlib import Path
from urllib.parse import urlparse
import os,requests
from requests.auth import HTTPBasicAuth

def download_urls(urls,output_dir,username=None,password=None,overwrite=False):
    """Download supplied direct file URLs with HTTP Range resume when supported.

    No credentials are stored. For PPS, use env PPS_USERNAME/PPS_PASSWORD.
    Public URLs need no authentication. Requires URLs to actual files.
    """
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    auth=HTTPBasicAuth(username,password) if username and password else None
    paths=[]
    with requests.Session() as session:
        for url in urls:
            if not url.startswith('https://'):raise ValueError('Only HTTPS URLs accepted')
            name=Path(urlparse(url).path).name
            if not name:raise ValueError('URL must point to a file')
            target=out/name;part=out/(name+'.part')
            if target.exists() and not overwrite:paths.append(str(target));continue
            offset=part.stat().st_size if part.exists() and not overwrite else 0
            headers={'Range':f'bytes={offset}-'} if offset else {}
            with session.get(url,auth=auth,headers=headers,stream=True,timeout=(30,120)) as response:
                response.raise_for_status()
                if offset and response.status_code!=206:offset=0
                mode='ab' if offset else 'wb'
                with part.open(mode) as f:
                    for chunk in response.iter_content(1024*1024):
                        if chunk:f.write(chunk)
            part.replace(target);paths.append(str(target))
    return paths

def download_pps(urls,output_dir,username=None,password=None):
    user=username or os.getenv('PPS_USERNAME');pw=password or os.getenv('PPS_PASSWORD')
    if not user:
        import getpass
        user=input('PPS registered username/email: ')
    if not pw:
        import getpass
        pw=getpass.getpass('PPS password: ')
    return download_urls(urls,output_dir,user,pw)
