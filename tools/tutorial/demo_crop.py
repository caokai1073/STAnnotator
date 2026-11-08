"""Crop the two SCP2927 sections to a window with mucosa, submucosa and muscle (cells + H&E + registration)."""
import sys, numpy as np, pandas as pd, anndata as ad, tifffile
# usage: demo_crop.py <folder with scp2927_<section>.h5ad, _HE.ome.tif, registration.csv> <out folder>  (needs anndata)
D=sys.argv[1].rstrip('/')+'/'
OUT=sys.argv[2]
WIN={'147328':(5000,8500,7000,11000),'167799':(1500,5500,0,5000)}
reg=pd.read_csv(D+'registration.csv')
rows=[]
for s,(x0,x1,y0,y1) in WIN.items():
    a=ad.read_h5ad(D+f'scp2927_{s}.h5ad')
    xy=a.obsm['spatial']; m=(xy[:,0]>=x0)&(xy[:,0]<x1)&(xy[:,1]>=y0)&(xy[:,1]<y1)
    b=a[m].copy()
    r=reg[reg.donor==int(s)].iloc[0]
    A=np.array([[r.a11,r.a12],[r.a21,r.a22]]); t=np.array([r.tx,r.ty])
    px=b.obsm['spatial']@A.T+t
    M=150  # margin in H&E pixels
    cx0=int(max(0,np.floor(px[:,0].min()-M))); cy0=int(max(0,np.floor(px[:,1].min()-M)))
    cx1=int(np.ceil(px[:,0].max()+M)); cy1=int(np.ceil(px[:,1].max()+M))
    import zarr
    store=tifffile.imread(D+f'scp2927_{s}_HE.ome.tif',aszarr=True,level=0); arr=zarr.open(store,mode='r')
    crop=np.asarray(arr[cy0:cy1,cx0:cx1,:]); store.close()
    print(s,'cells',b.n_obs,'H&E crop',crop.shape,'offset',cx0,cy0)
    with tifffile.TiffWriter(f'{OUT}/scp2927_{s}_HE_crop.ome.tif',bigtiff=False,ome=True) as tw:
        tw.write(crop,photometric='rgb',tile=(512,512),compression='jpeg',compressionargs={'level':88},subifds=4,
                 metadata={'axes':'YXS','PhysicalSizeX':0.69,'PhysicalSizeXUnit':'µm','PhysicalSizeY':0.69,'PhysicalSizeYUnit':'µm'})
        lv=crop
        for i in range(4):
            lv=lv[::2,::2]
            tw.write(np.ascontiguousarray(lv),photometric='rgb',tile=(512,512),compression='jpeg',compressionargs={'level':88},subfiletype=1)
    b.write_h5ad(f'{OUT}/scp2927_{s}_crop.h5ad',compression='gzip')
    rows.append(dict(slide=int(r.slide),donor=int(s),a11=r.a11,a12=r.a12,a21=r.a21,a22=r.a22,tx=r.tx-cx0,ty=r.ty-cy0))
pd.DataFrame(rows).to_csv(f'{OUT}/registration_crop.csv',index=False)
