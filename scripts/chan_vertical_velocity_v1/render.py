from pathlib import Path
import ast, json, hashlib, subprocess
import h5py, numpy as np, imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path.cwd()
# Read the frozen publication source map without executing its renderer.
tree=ast.parse((root/'scripts/chan_nov26_publication_v1/render.py').read_text())
cases=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='cases' for t in n.targets))
out=root/'rendered/chan-vertical-velocity-v1'
out.mkdir(parents=True,exist_ok=False)
cache={};limits=np.zeros(3)
for case,cfg in cases.items():
    times=list(range(0,cfg['end']+1,250));planes=[];sources=[]
    for t in times:
        part=next(root/'data'/p/'parent' for end,p in cfg['sources'] if t<=end)
        path=part/f'full_{t*1000:09d}.jld2'
        with h5py.File(path) as h:
            assert abs(float(h['t'][()])-t)<1e-7
            z=h['z'][:];indices=np.array([np.argmin(abs(z-a)) for a in (.2,.5,.9)])
            w=np.stack([.5*(h['w_face'][k,:,:]+h['w_face'][k+1,:,:]) for k in indices])
            assert w.shape==(3,256,256) and np.isfinite(w).all()
            # Check orientation and arithmetic centering against the saved diagnostic.
            if t==cfg['end']:
                wf=h['w_face'][:];wc=.5*(wf[:-1]+wf[1:])
                with h5py.File(part/f'frame_{t*1000:09d}.jld2') as frame:
                    np.testing.assert_allclose(np.sqrt(np.mean(wc.astype(float)**2,axis=(1,2))),frame['w_rms_profile'][:],rtol=2e-6,atol=1e-8)
                    np.testing.assert_allclose(z[indices],frame['video_heights'][:],atol=1e-12)
            planes.append(w);sources.append(str(path.relative_to(root)))
    data=np.array(planes);limits=np.maximum(limits,np.max(abs(data),axis=(0,2,3)))
    cache[case]=(cfg,times,data,z[indices],sources)
    print('READ_PASS',case,len(times),flush=True)
assert np.all(limits>0)
manifest={'status':'RENDERING','variable':'cell-centered vertical velocity w, positive upward','centering':'arithmetic average of the two adjacent native vertical faces','units':'native nondimensional velocity; no normalization by f or RMS','sampling_interval':250,'interpolation':'none; each saved snapshot held for one second','color_scale':'linear, symmetric; full sampled extrema per depth shared across all three cases and all times','limits_per_depth':limits.tolist(),'filtered':False,'cases':{}}
for case,(cfg,times,data,heights,sources) in cache.items():
    folder=out/case;folder.mkdir()
    fig,axs=plt.subplots(1,3,figsize=(16,5.3),dpi=100)
    fig.subplots_adjust(left=.05,right=.95,bottom=.14,top=.77,wspace=.5)
    images=[]
    for j,ax in enumerate(axs):
        im=ax.imshow(data[0,j],origin='lower',extent=[0,4,0,4],interpolation='nearest',cmap='RdBu_r',vmin=-limits[j],vmax=limits[j])
        images.append(im);ax.set(xlabel='x/H',ylabel='y/H',title=f'z/H = {heights[j]:.5f}')
        fig.colorbar(im,ax=ax,fraction=.045,pad=.03,label='w (upward positive)')
    title=fig.suptitle('',y=.97,fontsize=14)
    foot=fig.text(.5,.02,'Saved snapshots every Δt=250 · held 1 s each · no temporal interpolation · unfiltered · linear colors',ha='center',fontsize=10)
    movie=folder/f'vertical_velocity_0_{cfg["end"]}.mp4'
    writer=imageio_ffmpeg.write_frames(str(movie),(1600,530),fps=10,codec='libx264',bitrate='1200k',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=2,output_params=['-threads','2','-movflags','+faststart'])
    writer.send(None)
    stats=[]
    for i,t in enumerate(times):
        title.set_text(f"{cfg['title']} · Ω=3, β=0, no jet · t={t}\nVertical velocity · 256² × 80 · nondimensional time and velocity")
        for j,im in enumerate(images):im.set_data(data[i,j])
        fig.canvas.draw();rgb=np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()
        for repeat in range(10):writer.send(rgb)
        stats.append({'time':t,'min':data[i].min(axis=(1,2)).tolist(),'max':data[i].max(axis=(1,2)).tolist(),'rms':np.sqrt(np.mean(data[i].astype(float)**2,axis=(1,2))).tolist()})
    writer.close();fig.savefig(folder/'vertical_velocity_final.jpg',dpi=100,pil_kwargs={'quality':90});plt.close(fig)
    frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert frames==len(times)*10
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
    receipt={'case':case,'time_range':[0,cfg['end']],'saved_times':times,'heights':heights.tolist(),'source_files':sources,'centering_validation':'endpoint RMS profile agrees with saved simulation diagnostic','frames':frames,'fps':10,'duration_seconds':duration,'full_decode_passed':True,'snapshot_statistics':stats,'assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in folder.iterdir()}}
    manifest['cases'][case]=receipt
    print('MEDIA_PASS',case,frames,duration,flush=True)
manifest['status']='MEDIA_VERIFIED'
(out/'vertical_velocity_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('ALL_MEDIA_VERIFIED',flush=True)
