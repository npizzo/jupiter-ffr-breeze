"""Combine three depth videos into one synchronized row, without resampling time."""
import argparse,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
files=[a.folder/f'vorticity_{z}km.mp4' for z in (20,60,100)]
info=[]
for f in files:
 d=json.loads(subprocess.run(['ffprobe','-v','error','-select_streams','v:0','-show_entries','stream=width,height,nb_frames,r_frame_rate,duration','-of','json',str(f)],check=True,capture_output=True,text=True).stdout)['streams'][0];info.append(d)
if any(d!=info[0] for d in info[1:]):raise ValueError('Depth videos must have identical dimensions, frame count, cadence and duration')
a.output.parent.mkdir(parents=True,exist_ok=True)
cmd=['ffmpeg','-v','error','-y']
for f in files:cmd+=['-i',str(f)]
cmd+=['-filter_complex','[0:v][1:v][2:v]hstack=inputs=3[v]','-map','[v]','-c:v','libx264','-crf','24','-pix_fmt','yuv420p','-movflags','+faststart',str(a.output)]
subprocess.run(cmd,check=True)
subprocess.run(['ffmpeg','-v','error','-y','-sseof','-0.1','-i',str(a.output),'-frames:v','1',str(a.output.with_suffix('.png'))],check=True)
a.output.with_suffix('.json').write_text(json.dumps({'layout':'one row, three columns','depth_order_km':[20,60,100],'source_video_metadata':info,'temporal_resampling':False},indent=2)+'\n')
