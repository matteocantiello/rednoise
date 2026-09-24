# Regenerate all manuscript figures from saved artifacts. Run inside the astro kernel with:
#   exec(open("paper/make_figures.py").read())
import os, glob, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
try: apply_figure_style()
except NameError: pass
plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"legend.fontsize":6.5,"xtick.labelsize":7,"ytick.labelsize":7})
def _pl(ax,l):
    try: panel_letter(ax,l)
    except NameError: ax.text(-0.12,1.04,l,transform=ax.transAxes,fontweight="bold",fontsize=9)
OUT="paper/figures"; os.makedirs(OUT,exist_ok=True)
SC,DC=3.5,7.25
ell_sun=(5777.0)**4/(274*100)
spec_ell=lambda lT,lg: np.log10((10**np.asarray(lT))**4/(10**np.asarray(lg))/ell_sun)
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.pdf",bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.png",dpi=250,bbox_inches="tight")
    rr=fig.canvas.get_renderer()
    tx=[(t,t.get_window_extent(rr)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    tl=set()
    for a in fig.axes: tl|=set(a.get_xticklabels()+a.get_yticklabels())
    ov=[(a1.get_text()[:16],a2.get_text()[:16]) for i,(a1,b1) in enumerate(tx) for a2,b2 in tx[i+1:] if b1.overlaps(b2) and not(a1 in tl and a2 in tl)]
    print(f"{name}: overlaps={len(ov)} {ov[:3]}"); plt.close(fig)

# ---------------- data ----------------
EXT=pd.read_csv(host.artifact_path("114af09b-6e33-4161-a82e-8948b966202a"))
VM=pd.read_csv(host.artifact_path("dd29c64d-df65-4121-9718-725a60a3a01e"))
X=pd.read_parquet(host.artifact_path("000f7d1a-4ac6-4eec-81a1-3629b5611e2f"))
GR=np.load(host.artifact_path("b8244da2-cc76-42f7-8a0d-5fce81fd5677"))
ONS=pd.read_csv(host.artifact_path("ea99664e-a9aa-462f-bc51-f4d1a5da0aaa"))
EVO=pd.read_csv(host.artifact_path("bb4fef82-2184-4fc8-a891-23a8f0a854b0"))
ROT=pd.read_csv(host.artifact_path("9f5de9cc-6249-4398-a75d-4f4283c9cf6a"))
FIT=pd.read_csv(host.artifact_path("8bb974e0-2801-4ae8-92cb-b0f326559d01"))
P=X[(X["sample"]=="Bowman2020_Galactic")|((X["sample"]=="Shen2024_Galactic")&~X.dup_B2020_Shen)].copy()
P["log_vmac"]=np.log10(P.vmac.where(P.vmac>0)); P["log_vsini"]=np.log10(P.vsini.clip(lower=1))
EXT["log_alpha0"]=np.log10(EXT.alpha0.where(EXT.a0_ok&(EXT.alpha0>0))); EXT["log_nuchar"]=np.log10(EXT.nuchar)
ob=EXT[EXT.lT>=4.0].copy(); hot=VM[VM.logTeff_sp>4.3].copy(); cool=VM[VM.logTeff_sp<=4.3].copy()
# MIST tracks
trackdir="data/mist/MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS"
with open(f"{trackdir}/01500M.track.eep") as f:
    for l in f:
        if l.startswith("#") and "star_age" in l: hdr=l.strip("# \n").split(); break
ci={c:i for i,c in enumerate(hdr)}
tracks={}
for m in [7,9,12,15,20,24,32,40,60,85,120]:
    fn=f"{trackdir}/{int(m*100):05d}M.track.eep"
    if not os.path.exists(fn): continue
    a=np.genfromtxt(fn,comments="#"); k=(a[:,ci["phase"]]>=0)&(a[:,ci["phase"]]<=2)&(a[:,ci["log_Teff"]]>3.5)
    tracks[m]=dict(logT=a[k,ci["log_Teff"]],ell=spec_ell(a[k,ci["log_Teff"]],a[k,ci["log_g"]]))
def tracks_on(ax,lab=True,color="0.55",lw=0.6,xoff=0.006):
    for m,tk in tracks.items():
        ax.plot(tk["logT"],tk["ell"],color=color,lw=lw,zorder=1)
        if lab:
            i=np.argmin(tk["ell"]); ax.text(tk["logT"][i],tk["ell"][i]-0.05,f"{m}",fontsize=5,color=color,ha="center",va="top",zorder=1)
def bic_model(Xm,y):
    b,_,_,_=np.linalg.lstsq(Xm,y,rcond=None); r=y-Xm@b; n=len(y); return n*np.log((r**2).sum()/n)+Xm.shape[1]*np.log(n), b
def binned(x,y,edges):
    out=[]
    for a,b in zip(edges[:-1],edges[1:]):
        m=(x>=a)&(x<b)&np.isfinite(y)
        if m.sum()>=5: out.append((0.5*(a+b),np.median(y[m]),np.percentile(y[m],16),np.percentile(y[m],84),m.sum()))
    return pd.DataFrame(out,columns=["x","med","lo","hi","n"])
CS={"Bowman2020_Galactic":"#4053d3","Shen2024_Galactic":"#e6550d","Bowman2019b_LMC":"#2ca02c","Bowman2024_LMC":"#9467bd",
    "Bowman2024_SMC":"#8c564b","Ma2024_LMC_BSG":"#17becf","DornWallenstein2020_YSG":"#bcbd22","DornWallenstein2020_RSG":"#7f0000"}
LAB={"Bowman2020_Galactic":"Bowman+20 (MW)","Shen2024_Galactic":"Shen+24 (MW)","Bowman2019b_LMC":"Bowman+19b (LMC)","Bowman2024_LMC":"Bowman+24 (LMC)",
     "Bowman2024_SMC":"Bowman+24 (SMC)","Ma2024_LMC_BSG":"Ma+24 (LMC BSG)","DornWallenstein2020_YSG":"Dorn-Wallenstein+20 YSG","DornWallenstein2020_RSG":"Dorn-Wallenstein+20 RSG"}

# ---------------- Fig 1: samples on the sHRD ----------------
fig,axes=plt.subplots(1,2,figsize=(DC,3.6))
ax=axes[0]; tracks_on(ax)
for s,c in CS.items():
    d=EXT[EXT["sample"]==s]
    for tier,mk in [("A","o"),("B","s")]:
        dd=d[d.tier==tier]
        if len(dd): ax.scatter(dd.lT,dd.lL,s=14 if tier=="A" else 16,marker=mk,color=c,ec="k",lw=0.3,zorder=3,label=LAB[s] if tier=="A" or (d.tier=="A").sum()==0 else None)
ax.scatter([],[],marker="o",color="0.5",ec="k",lw=0.3,s=14,label="tier A (spectroscopic $\\log g$)"); ax.scatter([],[],marker="s",color="0.5",ec="k",lw=0.3,s=16,label="tier B (MIST-track mass)")
ax.set_xlim(4.80,3.45); ax.set_ylim(2.2,4.7); ax.set_xlabel(r"$\log T_{\rm eff}$ (K)"); ax.set_ylabel(r"$\log\,\mathcal{L}/\mathcal{L}_\odot$")
ax.legend(fontsize=5.2,loc="lower left",ncol=2,frameon=True,framealpha=0.9); ax.set_title(f"Red-noise samples (N={len(EXT)})")
ax=axes[1]; tracks_on(ax,lab=False)
sc=ax.scatter(VM.logTeff_sp,VM.logL_sp,c=VM.vmac,cmap="viridis",norm=mpl.colors.LogNorm(vmin=8,vmax=120),s=9,ec="k",lw=0.15,zorder=3)
cb=fig.colorbar(sc,ax=ax,pad=0.02); cb.set_label(r"$v_{\rm macro}$ (km s$^{-1}$)")
ax.set_xlim(4.72,4.0); ax.set_ylim(1.9,4.6); ax.set_xlabel(r"$\log T_{\rm eff}$ (K)"); ax.set_title(f"IACOB macroturbulence sample (N={len(VM)})")
for a,l in zip(axes,"ab"): _pl(a,l)
fig.tight_layout(); save(fig,"fig1_samples")

# ---------------- Fig 2: scalings on the primary sample ----------------
fig,axes=plt.subplots(2,2,figsize=(DC,5.8))
def panel(ax,x,y,xl,yl,logx=False,color_by=None,cl=None,legend_loc=None):
    sub=P[[x,y]+([color_by] if color_by else [])].dropna(); sub=sub[(sub[y]>0)&(sub[x]>0 if logx else True)]
    if color_by:
        sc=ax.scatter(sub[x],sub[y],c=sub[color_by],cmap="viridis",s=16,ec="k",lw=0.25,zorder=3); cb=fig.colorbar(sc,ax=ax,pad=0.02); cb.set_label(cl)
    else:
        for s,c,lab in [("Bowman2020_Galactic",CS["Bowman2020_Galactic"],"Bowman+20"),("Shen2024_Galactic",CS["Shen2024_Galactic"],"Shen+24")]:
            q=P[P["sample"]==s][[x,y]].dropna(); ax.scatter(q[x],q[y],s=16,color=c,ec="k",lw=0.25,zorder=3,label=lab)
        if legend_loc: ax.legend(loc=legend_loc)
    rho,p=spearmanr(sub[x],sub[y]); xx=np.log10(sub[x]) if logx else sub[x]; sl,ic=np.polyfit(xx,np.log10(sub[y]),1)
    xs=np.linspace(xx.min(),xx.max(),50); ax.plot(10**xs if logx else xs,10**(sl*xs+ic),color="0.2",lw=1.2,zorder=4)
    ax.set_yscale("log"); ax.set_xscale("log" if logx else "linear"); ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.text(0.03,0.96,f"$\\rho$={rho:+.2f}, p={p:.0e}\nN={len(sub)}, slope={sl:+.2f}",transform=ax.transAxes,va="top",fontsize=6.3,bbox=dict(fc="white",ec="none",alpha=0.8,pad=1.2))
panel(axes[0,0],"logTeff","nuchar",r"$\log T_{\rm eff}$ (K)",r"$\nu_{\rm char}$ (d$^{-1}$)",legend_loc="lower left"); axes[0,0].invert_xaxis()
panel(axes[0,1],"logL_spec","alpha0",r"$\log\,\mathcal{L}/\mathcal{L}_\odot$",r"$\alpha_0$ ($\mu$mag)",legend_loc="lower right")
panel(axes[1,0],"vmac","alpha0",r"$v_{\rm macro}$ (km s$^{-1}$)",r"$\alpha_0$ ($\mu$mag)",logx=True,color_by="logL_spec",cl=r"$\log\mathcal{L}$")
panel(axes[1,1],"tau","alpha0",r"$\tau$ (fractional MS age; Shen+24)",r"$\alpha_0$ ($\mu$mag)",color_by="logL_spec",cl=r"$\log\mathcal{L}$")
for a,l in zip(axes.ravel(),"abcd"): _pl(a,l)
fig.tight_layout(); save(fig,"fig2_scalings")

# ---------------- Fig 3: GP maps ----------------
fig,axes=plt.subplots(2,4,figsize=(DC,4.9))
SPEC=[("vmac",r"$v_{\rm macro}$ (km s$^{-1}$)","viridis",True,"OB"),("alpha0",r"$\alpha_0$ ($\mu$mag)","plasma",True,"OB"),
      ("nuchar",r"$\nu_{\rm char}$ (d$^{-1}$)","cividis",True,"ext"),("gamma",r"$\gamma$","magma",False,"ext")]
for j,(p,cl,cmap,logp,grid) in enumerate(SPEC):
    gx,gy=(GR["logTeff_OB"],GR["logL_OB"]) if grid=="OB" else (GR["logTeff_ext"],GR["logL_ext"])
    M,Sd,ok=GR[f"{p}_mean"],GR[f"{p}_std"],GR[f"{p}_ok"]
    ax=axes[0,j]; cf=ax.contourf(gx,gy,M,levels=12,cmap=cmap,zorder=1)
    ax.contourf(gx,gy,np.where(~ok,1.0,np.nan),levels=[0.5,1.5],colors=["0.92"],zorder=3)
    ax.contourf(gx,gy,np.where(~ok,1.0,np.nan),levels=[0.5,1.5],colors=["none"],hatches=["///"],zorder=3)
    cb=fig.colorbar(cf,ax=ax,location="top",orientation="horizontal",fraction=0.07,pad=0.03,aspect=22)
    if logp:
        t_=cb.get_ticks(); t_=t_[::2] if len(t_)>6 else t_; cb.set_ticks(t_); cb.set_ticklabels([f"{10**t:.0f}" if 10**t>=10 else f"{10**t:.2g}" for t in t_])
    cb.ax.tick_params(labelsize=5.5,pad=1); cb.set_label(cl,fontsize=6.5,labelpad=2)
    tracks_on(ax,lab=False,color="0.3")
    ax2=axes[1,j]; cf2=ax2.contourf(gx,gy,Sd,levels=10,cmap="Greys",zorder=1); ax2.contour(gx,gy,ok.astype(float),levels=[0.5],colors=["#d62728"],linewidths=1.0,zorder=3)
    cb2=fig.colorbar(cf2,ax=ax2,location="bottom",orientation="horizontal",fraction=0.07,pad=0.22,aspect=22); cb2.ax.tick_params(labelsize=5.5,pad=1); cb2.set_label("posterior s.d."+(" (dex)" if logp else ""),fontsize=6.5,labelpad=2)
    tracks_on(ax2,lab=False,color="0.3")
    xl=(4.72,4.03) if grid=="OB" else (4.78,3.50); yl=(1.9,4.5) if grid=="OB" else (1.9,4.7)
    for a in (ax,ax2): a.set_xlim(*xl); a.set_ylim(*yl); a.tick_params(labelsize=6)
    ax2.set_xlabel(r"$\log T_{\rm eff}$ (K)",labelpad=1)
    if grid=="ext":
        for a in (ax,ax2): a.axvline(4.0,color="w" if a is ax else "0.5",ls="--",lw=0.7,zorder=5)
for a in axes[:,0]: a.set_ylabel(r"$\log\,\mathcal{L}/\mathcal{L}_\odot$")
for a,l in zip(axes.ravel(),"abcdefgh"): _pl(a,l)
fig.subplots_adjust(left=0.06,right=0.99,top=0.93,bottom=0.12,wspace=0.28,hspace=0.25); save(fig,"fig3_gpmaps")

# ---------------- Fig 4: onset ----------------
fig,axes=plt.subplots(1,3,figsize=(DC,2.7))
FECZ=2.5
def panel_onset(ax,x,y,t,title,ylab,color,hinges,ylim,onset_txt):
    m=np.isfinite(x)&np.isfinite(y)&np.isfinite(t); x,y,t=x[m],y[m],t[m]
    ax.scatter(x,10**y,s=5,color=color,alpha=0.3,lw=0,zorder=2)
    b=binned(x,10**y,np.arange(2.0,4.6,0.2))
    ax.errorbar(b.x,b.med,yerr=[b.med-b.lo,b.hi-b.med],fmt="o",ms=3.5,color="k",ecolor="0.4",capsize=1.5,lw=0.8,zorder=4)
    xs=np.linspace(x.min(),x.max(),200); tm=np.median(t)
    if len(hinges)==1:
        _,cb=bic_model(np.column_stack([np.ones_like(x),x-hinges[0],np.clip(x-hinges[0],0,None),t]),y)
        ys=cb[0]+cb[1]*(xs-hinges[0])+cb[2]*np.clip(xs-hinges[0],0,None)+cb[3]*tm
    else:
        g1,g2=hinges; _,cb=bic_model(np.column_stack([np.ones_like(x),x-g1,np.clip(x-g1,0,None),np.clip(x-g2,0,None),t]),y)
        ys=cb[0]+cb[1]*(xs-g1)+cb[2]*np.clip(xs-g1,0,None)+cb[3]*np.clip(xs-g2,0,None)+cb[4]*tm
    ax.plot(xs,10**ys,color="#d62728",lw=1.4,zorder=5)
    for h in hinges: ax.axvline(h,color="#d62728",ls=":",lw=0.8)
    ax.axvline(FECZ,color="0.3",ls="--",lw=0.8)
    ax.set_yscale("log"); ax.set_ylim(*ylim); ax.set_xlim(2.0,4.5); ax.set_xlabel(r"$\log\,\mathcal{L}/\mathcal{L}_\odot$"); ax.set_ylabel(ylab); ax.set_title(title)
    ax.text(FECZ-0.03,ylim[0]*1.15,"FeCZ\nappears",fontsize=5.5,color="0.3",ha="right",va="bottom")
    ax.text(hinges[0]+0.03,ylim[0]*1.15,onset_txt,fontsize=5.5,color="#d62728",va="bottom")
panel_onset(axes[0],hot.logL_sp.values,np.log10(hot.vmac.values),hot.logTeff_sp.values,r"$v_{\rm macro}$, $\log T_{\rm eff}>4.3$ (N=579)",r"$v_{\rm macro}$ (km s$^{-1}$)","#4053d3",[3.16,3.70],(2,300),"onset 3.2 | sat. 3.7")
panel_onset(axes[1],cool.logL_sp.values,np.log10(cool.vmac.values),cool.logTeff_sp.values,r"$v_{\rm macro}$, $4.0<\log T_{\rm eff}\leq4.3$ (N=253)",r"$v_{\rm macro}$ (km s$^{-1}$)","#2ca02c",[2.78],(2,300),"hinge 2.8\n(n.s.)")
panel_onset(axes[2],ob.lL.values,ob.log_alpha0.values,ob.lT.values,r"$\alpha_0$, OB stars (N=222)",r"$\alpha_0$ ($\mu$mag)","#e6550d",[3.16],(4,3e4),"onset 3.2")
for a,l in zip(axes,"abc"): _pl(a,l)
fig.tight_layout(); save(fig,"fig4_onset")

# ---------------- Fig 5: evolution at fixed mass ----------------
obms=ob[(ob.ev_tau<1.0)&(ob.ev_dist<0.5)].copy(); VMms=VM[(VM.ev_tau<1.0)&(VM.ev_dist<0.5)].copy(); VMms["log_vmac"]=np.log10(VMms.vmac)
fig,axes=plt.subplots(1,3,figsize=(DC,2.7))
def cf_(t): r=EVO[(EVO.target==t)&(EVO.predictor=="ev_tau")].iloc[0]; return (r.coef,r.lo,r.hi)
def evo_panel(ax,df,y,ylab,title,coef,labels=True):
    dd=df.dropna(subset=[y,"ev_tau","ev_logM"])
    sc=ax.scatter(dd.ev_tau,10**dd[y],c=dd.ev_logM,cmap="plasma",s=9,ec="k",lw=0.2,zorder=3,vmin=0.8,vmax=1.8)
    m=LinearRegression().fit(dd[["ev_tau","ev_logM"]],dd[y]); tt=np.linspace(0,1,50)
    for lm in (1.0,1.3,1.6):
        ax.plot(tt,10**m.predict(np.column_stack([tt,np.full_like(tt,lm)])),color=plt.cm.plasma((lm-0.8)/1.0),lw=1.3,zorder=4)
        if labels: ax.text(1.02,10**m.predict([[1.0,lm]])[0],f"{10**lm:.0f}",fontsize=5.5,va="center",color=plt.cm.plasma((lm-0.8)/1.0))
    if not labels: ax.text(1.02,10**m.predict([[1.0,1.3]])[0],"10–40",fontsize=5.5,va="center",color="0.3")
    ax.set_yscale("log"); ax.set_xlim(-0.02,1.12); ax.set_xlabel(r"$\tau$ (fractional MS age)"); ax.set_ylabel(ylab); ax.set_title(title)
    ax.text(0.03,0.96,f"$\\partial\\log/\\partial\\tau$={coef[0]:+.2f} [{coef[1]:+.2f},{coef[2]:+.2f}]\nN={len(dd)}",transform=ax.transAxes,va="top",fontsize=6,bbox=dict(fc="white",ec="none",alpha=0.8,pad=1.2))
    return sc
sc=evo_panel(axes[0],obms,"log_alpha0",r"$\alpha_0$ ($\mu$mag)",r"$\alpha_0$",cf_("log_alpha0"))
evo_panel(axes[1],obms,"log_nuchar",r"$\nu_{\rm char}$ (d$^{-1}$)",r"$\nu_{\rm char}$",cf_("log_nuchar"),labels=False)
evo_panel(axes[2],VMms,"log_vmac",r"$v_{\rm macro}$ (km s$^{-1}$)",r"$v_{\rm macro}$",cf_("log_vmac"))
fig.subplots_adjust(left=0.07,right=0.905,top=0.90,bottom=0.17,wspace=0.45)
cax=fig.add_axes([0.925,0.17,0.012,0.73]); cb=fig.colorbar(sc,cax=cax); cb.set_label(r"$\log M/M_\odot$ (MIST)")
for a,l in zip(axes,"abc"): _pl(a,l)
save(fig,"fig5_evolution")

# ---------------- Fig 6: rotation ----------------
d=ob.dropna(subset=["f_rot","log_nuchar","lL","lT"]).copy(); d=d[d.vsini>0]; ratio=np.log10(d.nuchar/d.f_rot); near=d[d.rotmod_candidate]
fig,axes=plt.subplots(1,3,figsize=(DC,2.7))
ax=axes[0]; sc=ax.scatter(d.f_rot,d.nuchar,c=d.lT,cmap="viridis",s=12,ec="k",lw=0.25,zorder=3); cb=fig.colorbar(sc,ax=ax,pad=0.02); cb.set_label(r"$\log T_{\rm eff}$")
xx=np.logspace(-2.6,0.1,50)
for k,ls,lab in [(1,"-",r"$\nu_{\rm char}=f_{\rm rot}$"),(2,"--",r"$2f_{\rm rot}$"),(16.7,":",r"$16.7\,f_{\rm rot}$")]: ax.plot(xx,k*xx,color="0.3",ls=ls,lw=0.8,label=lab)
ax.scatter(near.f_rot,near.nuchar,s=40,fc="none",ec="#d62728",lw=0.9,zorder=4,label=f"near $f_{{\\rm rot}}$/$2f_{{\\rm rot}}$ (n={len(near)})")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(2e-3,1.5); ax.set_ylim(0.08,15); ax.legend(loc="upper left",fontsize=5.3)
ax.set_xlabel(r"$f_{\rm rot}=v\sin i/2\pi R$ (d$^{-1}$)"); ax.set_ylabel(r"$\nu_{\rm char}$ (d$^{-1}$)"); ax.set_title(f"N={len(d)}")
ax=axes[1]; ax.hist(ratio,bins=np.arange(-0.4,2.4,0.1),color="#4053d3",alpha=0.75,ec="k",lw=0.3)
for k,ls in [(0,"-"),(np.log10(2),"--")]: ax.axvline(k,color="0.3",ls=ls,lw=0.8)
ax.axvline(np.median(ratio),color="#d62728",lw=1.2); ax.text(np.median(ratio)+0.05,ax.get_ylim()[1]*0.9,f"median {10**np.median(ratio):.0f}×",color="#d62728",fontsize=6)
ax.set_xlabel(r"$\log_{10}(\nu_{\rm char}/f_{\rm rot})$"); ax.set_ylabel("stars"); ax.set_title("no pile-up at $f_{\\rm rot}$, $2f_{\\rm rot}$")
ax=axes[2]; rows=[(r"$\nu_{\rm char}$",ROT.loc[ROT.quantity=="log_nuchar"].iloc[0]),(r"$\alpha_0$",ROT.loc[ROT.quantity=="log_alpha0"].iloc[0]),(r"$v_{\rm macro}$",ROT.loc[ROT.quantity=="log_vmac"].iloc[0])]
ys=np.arange(3)[::-1]
for y_,(lab,r) in zip(ys,rows):
    ax.plot([0,r.slope_log_vsini],[y_,y_],color="0.6",lw=2,zorder=2); ax.scatter(r.slope_log_vsini,y_,s=50,color="#4053d3",ec="k",lw=0.5,zorder=3)
    ax.text(r.slope_log_vsini+0.04,y_,f"$\\rho$={r.partial_rho:+.2f}, p={r.p:.0e}\nN={int(r.N)}",va="center",fontsize=5.6)
ax.axvline(0,color="k",lw=0.7); ax.axvline(1,color="0.3",ls="--",lw=0.8); ax.text(1.0,ys[0]+0.45,"∝ rotation rate",ha="center",fontsize=5.5,color="0.3")
ax.set_yticks(ys); ax.set_yticklabels([r_[0] for r_ in rows]); ax.set_xlim(-0.1,1.3); ax.set_ylim(-0.6,2.9)
ax.set_xlabel(r"$\partial\log(\cdot)/\partial\log v\sin i$ after $\mathcal{L},T_{\rm eff}$"); ax.set_title("weak common multiplier")
for a,l in zip(axes,"abc"): _pl(a,l)
fig.tight_layout(); save(fig,"fig6_rotation")

# ---------------- Fig 7: theory confrontation ----------------
fig=plt.figure(figsize=(DC,2.9)); gs=fig.add_gridspec(1,3,width_ratios=[1.35,0.8,1.0],wspace=0.45)
axa=fig.add_subplot(gs[0,0])
tests=[(r"$\partial\nu_{\rm char}/\partial T_{\rm eff}$","+1.53",+1,+1,-1),
       (r"$\partial\alpha_0/\partial\tau$ at fixed $M$","+1.73 (×50)",+1,+1,-1),
       (r"$\alpha_0$–$v_{\rm macro}$ (partial)","+0.28",+1,+1,-1),
       (r"onset at $\log\mathcal{L}\approx3.0$–3.2","flat below",+1,+1,-1),
       (r"amplitude budget","350 vs 0.1 $\mu$mag",+1,+1,-1),
       (r"$\nu_{\rm char}$ vs $f_{\rm rot}$","17×, weak",0,0,0),
       (r"LMC vs SMC $\nu_{\rm char}$","equal (method-limited)",0,0,0)]
for i,(nm,obs,osign,fs,isg) in enumerate(tests):
    y=len(tests)-i
    axa.text(-0.02,y,nm,ha="right",va="center",fontsize=6.2); axa.text(0.27,y,obs,ha="center",va="center",fontsize=5.8)
    for x,sgn in [(0.62,fs),(0.87,isg)]:
        if osign==0: col,mark="0.75","○"
        else: match=(sgn==osign); col,mark=("#2ca25f" if match else "#de2d26"),("✓" if match else "✗")
        axa.scatter(x,y,s=95,color=col,ec="k",lw=0.4,zorder=3); axa.text(x,y,mark,ha="center",va="center",fontsize=6.5,color="white" if osign else "0.3",zorder=4)
axa.text(0.27,len(tests)+1,"observed",ha="center",fontsize=6.5,fontweight="bold"); axa.text(0.62,len(tests)+1,"FeCZ",ha="center",fontsize=6.5,fontweight="bold"); axa.text(0.87,len(tests)+1,"core IGW",ha="center",fontsize=6.5,fontweight="bold")
axa.set_xlim(-0.05,1.0); axa.set_ylim(0.3,len(tests)+1.6); axa.axis("off"); axa.set_title("Diagnostic scorecard",fontsize=8.5)
axb=fig.add_subplot(gs[0,1])
a0=P.alpha0.dropna(); lo,med,hi=a0.quantile([.1,.5,.9])
axb.plot([0,0],[lo,hi],color="#4053d3",lw=3,solid_capstyle="round"); axb.scatter(0,med,s=60,color="#4053d3",ec="k",zorder=3)
axb.plot([1,1],[0.06,0.21],color="#de2d26",lw=3,solid_capstyle="round"); axb.scatter(1,0.06,s=40,color="#de2d26",ec="k",zorder=3)
axb.set_yscale("log"); axb.set_ylim(0.02,3e4); axb.set_xlim(-0.6,1.6); axb.set_xticks([0,1]); axb.set_xticklabels(["observed\n(10–90%)","core IGW\n(Anders+23)"],fontsize=6)
axb.set_ylabel(r"$\alpha_0$ ($\mu$mag)"); axb.set_title("Amplitude budget",fontsize=8.5)
axb.annotate("",xy=(0.5,0.15),xytext=(0.5,med),arrowprops=dict(arrowstyle="<->",color="0.3",lw=0.8)); axb.text(0.55,8,f"×{med/0.15:.0f}",fontsize=6.5,color="0.3")
axc=fig.add_subplot(gs[0,2])
dd=P[["log_alpha0","logL_spec","logTeff","log_vmac"]].dropna()
ry=dd.log_alpha0-LinearRegression().fit(dd[["logL_spec","logTeff"]],dd.log_alpha0).predict(dd[["logL_spec","logTeff"]])
rx=dd.log_vmac-LinearRegression().fit(dd[["logL_spec","logTeff"]],dd.log_vmac).predict(dd[["logL_spec","logTeff"]])
rp,pp=spearmanr(rx,ry); sl=np.polyfit(rx,ry,1)[0]
axc.scatter(rx,ry,s=12,color="#4053d3",ec="k",lw=0.25,alpha=0.8); xs=np.linspace(rx.min(),rx.max(),20); axc.plot(xs,np.polyval(np.polyfit(rx,ry,1),xs),color="#d62728",lw=1.2)
axc.axhline(0,color="grey",ls=":",lw=0.7); axc.axvline(0,color="grey",ls=":",lw=0.7)
axc.set_xlabel(r"$\log v_{\rm macro}$ residual"); axc.set_ylabel(r"$\log\alpha_0$ residual"); axc.set_title("Added-variable test",fontsize=8.5)
axc.text(0.03,0.96,f"$\\rho$={rp:+.2f}, p={pp:.3f}\nN={len(dd)}, slope={sl:+.2f}",transform=axc.transAxes,va="top",fontsize=6,bbox=dict(fc="white",ec="none",alpha=0.8,pad=1.2))
for a,l in zip([axa,axb,axc],"abc"): _pl(a,l)
save(fig,"fig7_theory")
print("added-variable:",round(rp,3),round(pp,4),len(dd),round(sl,3))

# ---------------- Appendix A1: Bowman vs Shen overlap ----------------
import re
def norm(n):
    s=str(n).upper().replace("−","-"); s=re.sub(r"\s+","",s); s=s.replace("HDE","HD"); return re.sub(r"^(HD|BD|CPD|ALS|HR)0+(\d)",r"\1\2",s)
b=X[X["sample"]=="Bowman2020_Galactic"].assign(key=lambda t:t.star.map(norm)); s=X[X["sample"]=="Shen2024_Galactic"].assign(key=lambda t:t.star.map(norm))
cmp_=b.merge(s,on="key",suffixes=("_B","_S"))
fig,axes=plt.subplots(1,3,figsize=(DC,2.5))
for ax,(c,lab,lg) in zip(axes,[("alpha0",r"$\alpha_0$ ($\mu$mag)",True),("nuchar",r"$\nu_{\rm char}$ (d$^{-1}$)",True),("gamma",r"$\gamma$",False)]):
    x,y=cmp_[c+"_B"],cmp_[c+"_S"]; ax.scatter(x,y,s=16,color="#4053d3",ec="k",lw=0.3,zorder=3)
    lim=[min(x.min(),y.min())*0.8,max(x.max(),y.max())*1.25] if lg else [min(x.min(),y.min())-0.3,max(x.max(),y.max())+0.3]
    ax.plot(lim,lim,"k--",lw=0.8); ax.set_xlim(lim); ax.set_ylim(lim)
    if lg: ax.set_xscale("log"); ax.set_yscale("log"); dlt=np.log10(y/x); unit="dex"
    else: dlt=y-x; unit=""
    rho,p=spearmanr(x,y)
    ax.text(0.03,0.96,f"N={len(cmp_)}\n$\\Delta$ median={dlt.median():+.2f} {unit}\nMAD={(dlt-dlt.median()).abs().median():.2f}\n$\\rho$={rho:.2f}",transform=ax.transAxes,va="top",fontsize=6,bbox=dict(fc="white",ec="none",alpha=0.8,pad=1.2))
    ax.set_xlabel(lab+" Bowman+20"); ax.set_ylabel(lab+" Shen+24")
for a,l in zip(axes,"abc"): _pl(a,l)
fig.tight_layout(); save(fig,"figA1_overlap")
print("overlap N:",len(cmp_))

# ---------------- Appendix A2: tier-B calibration ----------------
pts=[];lm=[]
for fn in sorted(glob.glob(f"{trackdir}/*M.track.eep")):
    m0=int(os.path.basename(fn)[:5])/100
    if m0<5 or m0>120: continue
    a=np.genfromtxt(fn,comments="#"); k=(a[:,ci["phase"]]>=0)&(a[:,ci["phase"]]<=6)
    pts.append(np.column_stack([a[k,ci["log_Teff"]],a[k,ci["log_L"]]])); lm.append(np.log10(a[k,ci["star_mass"]]))
PTS=np.vstack(pts); LM=np.concatenate(lm); SCL=np.array([0.02,0.10]); tree=cKDTree(PTS/SCL)
def mass_from_track(lT,lL):
    dst,i=tree.query(np.column_stack([lT,lL])/SCL,k=8); w=1/(dst+0.3)**2; return (w*LM[i]).sum(1)/w.sum(1)
magA=EXT[(EXT.tier=="A")&EXT["sample"].isin(["Bowman2024_LMC","Bowman2024_SMC","Ma2024_LMC_BSG"])].merge(X[["star","logL_class"]],on="star",how="left").dropna(subset=["logL_class"])
ell_est=magA.logL_class.values-mass_from_track(magA.lT.values,magA.logL_class.values)
dlt=ell_est-magA.lL.values
fig,ax=plt.subplots(figsize=(SC,2.9))
ax.scatter(magA.lL,ell_est,s=18,color="#9467bd",ec="k",lw=0.3,zorder=3); lim=[3.2,4.7]; ax.plot(lim,lim,"k--",lw=0.8); ax.plot(lim,np.array(lim)+np.median(dlt),color="#d62728",lw=1,label=f"offset {np.median(dlt):+.2f} dex")
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel(r"$\log\mathcal{L}$ from spectroscopic $\log g$ (tier A)"); ax.set_ylabel(r"$\log\mathcal{L}$ from $\log L-\log M_{\rm MIST}$ (tier-B method)")
ax.text(0.03,0.96,f"N={len(magA)}\nmedian {np.median(dlt):+.3f}, MAD {np.median(np.abs(dlt-np.median(dlt))):.3f}, s.d. {dlt.std():.3f} dex",transform=ax.transAxes,va="top",fontsize=6,bbox=dict(fc="white",ec="none",alpha=0.8,pad=1.2)); ax.legend(loc="lower right")
fig.tight_layout(); save(fig,"figA2_tierB")
print("tierB calib:",len(magA),round(np.median(dlt),3),round(dlt.std(),3))

# ---------------- Appendix B1: correlation matrix ----------------
rn=["log_nuchar","log_alpha0","gamma"]; sp=["logL_spec","logTeff","logg","vsini","vmac","tau","mass"]
rho=np.full((3,7),np.nan); pv=rho.copy(); NN=np.zeros_like(rho,dtype=int)
for i,a_ in enumerate(rn):
    for j,b_ in enumerate(sp):
        s_=P[[a_,b_]].dropna()
        if len(s_)>=10: rho[i,j],pv[i,j]=spearmanr(s_[a_],s_[b_]); NN[i,j]=len(s_)
fig,ax=plt.subplots(figsize=(DC*0.75,2.4)); im=ax.imshow(rho,cmap="RdBu_r",vmin=-1,vmax=1,aspect="auto")
ax.set_xticks(range(7)); ax.set_xticklabels([r"$\log\mathcal{L}$",r"$\log T_{\rm eff}$",r"$\log g$",r"$v\sin i$",r"$v_{\rm macro}$",r"$\tau$",r"$M$"])
ax.set_yticks(range(3)); ax.set_yticklabels([r"$\log\nu_{\rm char}$",r"$\log\alpha_0$",r"$\gamma$"])
for i in range(3):
    for j in range(7):
        if np.isnan(rho[i,j]): continue
        st="***" if pv[i,j]<1e-3 else "**" if pv[i,j]<1e-2 else "*" if pv[i,j]<0.05 else ""
        ax.text(j,i,f"{rho[i,j]:+.2f}{st}\nn={NN[i,j]}",ha="center",va="center",fontsize=5.8,color="white" if abs(rho[i,j])>0.55 else "black")
cb=fig.colorbar(im,ax=ax,pad=0.02); cb.set_label(r"Spearman $\rho$"); fig.tight_layout(); save(fig,"figB1_corr")
print("ALL FIGURES DONE")
