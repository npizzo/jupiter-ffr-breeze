"""Build fifty-day gallery assets on Unity, after both segment-five jobs succeed."""
from pathlib import Path
import tomllib,json,shutil,subprocess,os,hashlib
import imageio_ffmpeg
from PIL import Image
root=Path.cwd();source=root/'data/cascade100-v1'
work=root/'rendered/cascade100-first50-work';out=root/'rendered/cascade100-first50-gallery'
work.mkdir(exist_ok=True);out.mkdir(exist_ok=True)
provenance={'status':'PREPARING','end_day':50,'campaign_target_day':100,'expected_frames_per_movie':401,'compute_location':'Unity','raw_fields_downloaded':False,'cases':{},'renderer_sha256':hashlib.sha256((root/'scripts/render_cascade50_movies.py').read_bytes()).hexdigest()}
for mode in ['heated','control']:
    case=source/('cascade100_f_'+mode);view=work/mode;view.mkdir(exist_ok=True)
    records=[]
    for seg in range(1,6):
        original=case/f'segment_{seg}';m=tomllib.loads((original/'metadata.toml').read_text())
        assert m['status']=='COMPLETE' and m['stop_seconds']==seg*10*86400
        assert m['nx']==256 and m['ny']==256 and m['nz']==64
        assert m['shear_rossby']==0 and m['beta']==0
        assert json.loads((original/'scale_transfer/provenance.json').read_text())['status']=='COMPLETE'
        link=view/f'segment_{seg}'
        if not link.exists():link.symlink_to(original,target_is_directory=True)
        assert link.resolve()==original.resolve()
        shutil.copyfile(original/'metadata.toml',out/f'{mode}_segment{seg}_metadata.toml')
        records.append({'segment':seg,'source':str(original.relative_to(root)),'frames':len(list(original.glob('frame_*.jld2'))),'driver_sha256':m['script_sha256'],'configuration_sha256':m['configuration_sha256']})
    assert sum(item['frames'] for item in records)==401
    provenance['cases'][mode]=records
base=[str(root/'tools/gallery-env/bin/python'),str(root/'scripts/render_cascade50_movies.py'),str(work/'heated'),str(work/'control'),'--labels','Heating ON','Heating OFF','--title','Coarse f-plane | 1200-km seeds | no imposed shear | days 0–50','--tracer','--fps','12','--crf','29']
media={}
for name,index in [('comparison',None),('heated',0),('unheated',1)]:
    target=out/f'{name}_0_50days.mp4'
    cmd=base+['--output',str(target)]
    if index is not None:cmd+=['--case-index',str(index)]
    subprocess.run(cmd,check=True)
    manifest=json.loads(target.with_suffix('.json').read_text())
    assert manifest['frames']==401 and manifest['times_seconds'][0]==0 and manifest['times_seconds'][-1]==50*86400
    frames,duration=imageio_ffmpeg.count_frames_and_secs(str(target))
    assert frames==401 and abs(duration-401/12)<.03
    # GitHub browser upload accepts at most 25 MB/file. Re-encode on Unity only
    # if needed, preserving dimensions, frame rate, color interpretation and all frames.
    crf=31
    while target.stat().st_size>24_000_000 and crf<=39:
        temp=target.with_name(target.stem+'_web.mp4')
        subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(target),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        temp.replace(target);manifest['web_reencode_crf']=crf;crf+=2
    assert target.stat().st_size<24_000_000,'Movie too large for browser upload'
    frames,duration=imageio_ffmpeg.count_frames_and_secs(str(target));assert frames==401
    manifest['final_mp4_sha256']=hashlib.sha256(target.read_bytes()).hexdigest()
    target.with_suffix('.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with Image.open(target.with_suffix('.png')) as im:dimensions=list(im.size);im.verify()
    media[target.name]={'frames':frames,'duration_seconds':duration,'bytes':target.stat().st_size,'sha256':manifest['final_mp4_sha256'],'poster_pixels':dimensions}
    print('Verified',target.name,media[target.name],flush=True)
provenance['status']='COMPLETE';provenance['media']=media
(out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
print('READY_FOR_PUBLICATION',out,flush=True)
