import json, numpy as np, pandas as pd, anndata as ad
from scipy import ndimage
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
D=sys.argv[1].rstrip('/')+'/'  # usage: demo_geojson.py <out folder of demo_crop.py>: QuPath-style outlines
reg=pd.read_csv(D+'registration_crop.csv')
for s in ['147328','167799']:
    a=ad.read_h5ad(D+f'scp2927_{s}_crop.h5ad',backed='r')
    r=reg[reg.donor==int(s)].iloc[0]
    xy=a.obsm['spatial'][:]
    px=np.c_[xy[:,0]*r.a11+xy[:,1]*r.a12+r.tx, xy[:,0]*r.a21+xy[:,1]*r.a22+r.ty]
    feats=[]
    for lab,name,color in [('Submucosa','Submucosa',[46,139,87]),('Musc.Propria','Muscularis propria',[70,90,200])]:
        m=(a.obs['histopathology']==lab).values
        g=40; W=int(px[:,0].max()/g)+3; H=int(px[:,1].max()/g)+3
        grid=np.zeros((H,W)); np.add.at(grid,(np.clip((px[m,1]/g).astype(int),0,H-1),np.clip((px[m,0]/g).astype(int),0,W-1)),1)
        grid=ndimage.gaussian_filter(grid,2.0)
        cs=plt.contour(grid,levels=[0.35])
        polys=[]
        for path in cs.get_paths():
            for poly in path.to_polygons():
                if len(poly)<20: continue
                area=0.5*abs(np.dot(poly[:,0],np.roll(poly[:,1],1))-np.dot(poly[:,1],np.roll(poly[:,0],1)))
                if area<400: continue
                pts=(poly*g).round(1).tolist()
                if pts[0]!=pts[-1]: pts.append(pts[0])
                polys.append(pts[::3]+[pts[0]])
        plt.close('all')
        polys.sort(key=len, reverse=True)
        for p in polys[:3]:
            feats.append({"type":"Feature","id":f"{name}-{len(feats)}","geometry":{"type":"Polygon","coordinates":[p]},
                          "properties":{"objectType":"annotation","classification":{"name":name,"color":color}}})
        print(s,name,len(polys),[len(p) for p in polys[:3]])
    json.dump({"type":"FeatureCollection","features":feats},open(D+f'scp2927_{s}_qupath.geojson','w'))
