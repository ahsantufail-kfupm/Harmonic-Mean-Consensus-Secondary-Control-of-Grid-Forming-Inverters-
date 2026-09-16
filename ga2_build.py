import os, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

# ---------------------------------------------------------------------------
# Frequency traces for the right-hand panel: regenerated on the fly (numpy
# only, ~2 s) with the same 6-unit test system, RK4 @ 1 ms, Ts = 20 ms and
# conditional anti-windup as paper1_reproduce.m. Feasible: +30 MW at bus 4;
# infeasible: +75 MW at bus 3 (0.27 pu deficit -> -0.233 Hz/s, Prop. 4).
# ---------------------------------------------------------------------------
def _ga_traces():
    f0=50.0; w0=2*np.pi*f0; N=6
    M=2*np.array([5,4,6,4.5,5.5,4.])/w0; Dv=1.2*np.ones(N)
    Pmax=np.array([.30,.25,.20,.20,.25,.20]); Pmin=np.zeros(N)
    KDR=np.array([.20,.18,.20,.18,.20,.18])
    Pl0=np.array([.08,.16,.12,.12,.14,.10])
    AL,KS,TS=0.30,0.8,0.02; RX=0.3; X=1/8; R=RX*X; G_L=R/(R*R+X*X); B_L=X/(R*R+X*X)
    LINES=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(0,2),(1,4)]
    COMM =[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(0,2),(0,3)]
    gen=np.array([1,1,1,1,1,0],bool)
    dg=np.zeros(N,int)
    for i,j in COMM:
        if gen[i] and gen[j]: dg[i]+=1; dg[j]+=1
    A=np.zeros((N,N))
    for i,j in COMM:
        if gen[i] and gen[j]:
            w=(dg[i]+dg[j])/(2*dg[i]*dg[j]); A[i,j]=w; A[j,i]=w
    ep=0.6/A[gen].sum(1).max()
    W=np.eye(N); W[np.ix_(gen,gen)]=(np.eye(gen.sum())
        - ep*(np.diag(A[np.ix_(gen,gen)].sum(1))-A[np.ix_(gen,gen)]))
    def run(T,bus,dP):
        dt=1e-3; steps=int(T/dt); ns=int(TS/dt)
        Pl=Pl0.copy(); de=np.zeros(N); dw=np.zeros(N)
        r0=Pl[gen].sum()/Pmax[gen].sum(); u=r0*Pmax*gen
        inj=lambda uu,dd: np.clip(uu-KDR*dd,Pmin,Pmax)*gen
        F=np.zeros(steps); t=np.arange(steps)*dt
        def f(dE,dW):
            Pi=inj(u,dW); pe=np.zeros(N)
            for i,j in LINES:
                th=dE[i]-dE[j]
                pe[i]+=G_L*(1-np.cos(th))+B_L*np.sin(th)
                pe[j]+=G_L*(1-np.cos(th))-B_L*np.sin(th)
            damp=Dv*(dW-dW.mean())*gen
            return dW.copy(),(Pi-Pl-pe-damp)/M
        for k in range(steps):
            if k*dt>=1.0 and Pl[bus]==Pl0[bus]: Pl[bus]+=dP
            a1,b1=f(de,dw); a2,b2=f(de+dt/2*a1,dw+dt/2*b1)
            a3,b3=f(de+dt/2*a2,dw+dt/2*b2); a4,b4=f(de+dt*a3,dw+dt*b3)
            de+=dt/6*(a1+2*a2+2*a3+a4); dw+=dt/6*(b1+2*b2+2*b3+b4)
            if k%ns==0:
                Pi=inj(u,dw); r=np.zeros(N); r[gen]=Pi[gen]/Pmax[gen]
                rm=(W@r); du=AL*(-dw)*TS+KS*(rm-r)*Pmax; du*=gen
                hi=(Pi>=Pmax-1e-6)&(du>0); lo=(Pi<=Pmin+1e-6)&(du<0)
                du[hi]=0; du[lo]=0; u=np.clip(u+du,Pmin,Pmax)
            F[k]=f0+dw[gen].mean()/(2*np.pi)
        return t,F
    t,fF=run(8.0,3,0.30); _,fI=run(8.0,2,0.75)
    return dict(t=t,fF=fF,fI=fI)

if os.path.exists('ga_data.pkl'):
    d=pickle.load(open('ga_data.pkl','rb'))
else:
    d=_ga_traces()

BLUE='#1a5fa8'; LBLUE='#5b9bd5'; DGRAY='#555555'; RED='#c0392b'; GREEN='#1e7d32'
GOLD='#b8860b'; CARD='#f5f8fc'; BG='white'
plt.rcParams.update({'font.family':['Liberation Sans','DejaVu Sans'],'font.size':12,
                     'axes.edgecolor':'#999999','axes.linewidth':0.9})
fig=plt.figure(figsize=(13.28,5.31),dpi=300); fig.patch.set_facecolor(BG)

# ================= header =================
fig.text(0.5,0.965,'Harmonic-Mean Consensus Secondary Control of Grid-Forming Inverters',
         ha='center',va='center',fontsize=16.5,fontweight='bold',color='#1c1c1c')
fig.text(0.5,0.915,'each inverter shares one number with its neighbours  '
         r'$\rightarrow$  provably overload-free power sharing and a closed-form recoverability clock',
         ha='center',va='center',fontsize=11.8,color=DGRAY,style='italic')

# big flow arrows
for x0 in (0.345,0.612):
    fig.patches.append(FancyArrowPatch((x0,0.50),(x0+0.020,0.50),transform=fig.transFigure,
        arrowstyle='-|>',mutation_scale=34,lw=3,color=BLUE))

# ================= Panel A: system + law =================
axA=fig.add_axes([0.015,0.115,0.325,0.755]); axA.set_xlim(0,10); axA.set_ylim(0,10); axA.axis('off')
axA.text(5,9.55,'ISLANDED AC MICROGRID',ha='center',fontsize=12,fontweight='bold',color=BLUE)
pos={1:(5,7.9),2:(8.2,5.9),3:(7.1,2.6),4:(2.9,2.6),5:(1.8,5.9)}
elec=[(1,2),(2,3),(3,4),(4,5),(1,3),(2,5)]                 # electrical G_E (paper Fig. 1)
comm=[(1,2),(2,3),(3,4),(4,5),(1,3),(1,4)]                 # communication G_C (distinct)
for i,j in elec:
    axA.plot(*zip(pos[i],pos[j]),color='#9aa0a6',lw=2.6,zorder=1,solid_capstyle='round')
for i,j in comm:
    axA.add_patch(FancyArrowPatch(pos[i],pos[j],connectionstyle='arc3,rad=0.22',
        arrowstyle='-',ls=(0,(4,3)),lw=1.8,color=LBLUE,zorder=2))
# exchange arrows on one highlighted comm arc (1-2)
axA.add_patch(FancyArrowPatch(pos[1],pos[2],connectionstyle='arc3,rad=0.22',
    arrowstyle='<|-|>',mutation_scale=14,lw=1.8,color=BLUE,zorder=3))
mid=(6.95,7.45); axA.text(*mid,r'$r_j$',fontsize=12,color=BLUE,fontweight='bold')
for i,(x,y) in pos.items():
    axA.add_patch(FancyBboxPatch((x-0.72,y-0.52),1.44,1.04,boxstyle='round,pad=0.06',
        fc='white',ec=BLUE,lw=2.2,zorder=4))
    axA.plot([x-0.42,x-0.12],[y+0.13,y+0.13],color=BLUE,lw=1.6,zorder=5)          # DC bar
    th=np.linspace(0,2*np.pi,60)
    axA.plot(x+0.24+0.20*th/6.28*0+0.20*np.linspace(-1,1,60),
             y+0.13+0.14*np.sin(np.linspace(0,2*np.pi,60)),color=BLUE,lw=1.4,zorder=5)  # AC sine
    axA.text(x,y-0.27,rf'$r_{i}$',ha='center',fontsize=11.5,color='#1c1c1c',zorder=5)
# load arrow at bus 3
axA.add_patch(FancyArrowPatch((pos[3][0]+0.2,pos[3][1]-0.60),(pos[3][0]+0.2,pos[3][1]-1.25),
    arrowstyle='-|>',mutation_scale=15,lw=2.2,color=DGRAY,zorder=3))
axA.text(pos[3][0]+0.52,pos[3][1]-0.95,'load',fontsize=10,color=DGRAY,va='center')
# legend
axA.plot([0.35,1.05],[9.32,9.32],color='#9aa0a6',lw=2.6); axA.text(1.22,9.32,'electrical network',fontsize=9.6,va='center',color=DGRAY)
axA.plot([0.35,1.05],[8.80,8.80],color=LBLUE,lw=1.8,ls=(0,(4,3))); axA.text(1.22,8.80,'communication (harmonic weights)',fontsize=9.6,va='center',color=DGRAY)
# law box
axA.text(5,0.18,r'$\dot u_i=\alpha(\omega^{\star}\!-\omega_i)+\beta\sum_{j\in N_i} a_{ij}(r_j-r_i)$'
                 '\n'+r'$a_{ij}=(d_i\!+\!d_j)\,/\,(2d_id_j),\qquad r_i=P_i\,/\,P_i^{\max}$',
         ha='center',va='center',fontsize=11.4,color='#1c1c1c',linespacing=1.5,
         bbox=dict(boxstyle='round,pad=0.35',fc=CARD,ec=BLUE,lw=1.4))

# ================= Panel B: three guarantee cards =================
axB=fig.add_axes([0.372,0.075,0.235,0.80]); axB.set_xlim(0,10); axB.set_ylim(0,10); axB.axis('off')
cards=[(8.35,GOLD ,'THEOREM 1 — fastest classical weighting',
        r'$\lambda_2^{\mathrm{harm}}\geq\lambda_2^{\mathrm{geom}}\geq\lambda_2^{\mathrm{arith}}\geq\lambda_2^{\mathrm{quad}}$',
        'proven ordering, exact equality condition'),
       (5.15,GREEN,'THEOREM 3 — minimax-safe sharing',
        r'$r_i=r^{\star}=\sum_i P_i^{\ell}\,/\sum_i P_i^{\max}<1$',
        'no unit pinned; peak inverter stress minimal'),
       (1.95,RED  ,'PROPOSITION 4 — recoverability limit',
        r'$\dot\omega_M=-\Delta P\,/\sum_i M_i$  once infeasible',
        'seconds to load shedding known in advance')]
for yc,c,t1,t2,t3 in cards:
    axB.add_patch(FancyBboxPatch((0.12,yc-1.42),9.76,2.84,boxstyle='round,pad=0.12',
        fc=CARD,ec='#c8d4e2',lw=1.2))
    axB.add_patch(Rectangle((0.12,yc-1.42),0.32,2.84,fc=c,ec='none'))
    axB.text(0.78,yc+0.84,t1,fontsize=9.6,fontweight='bold',color='#1c1c1c',va='center')
    axB.text(0.78,yc-0.05,t2,fontsize=11.4,color=c,va='center')
    axB.text(0.78,yc-0.92,t3,fontsize=9.8,color=DGRAY,va='center')

# ================= Panel C: evidence (real data) =================
axC1=fig.add_axes([0.665,0.545,0.315,0.315])
units=np.arange(5); Pmax=np.array([.30,.25,.20,.20,.25]); r_mw=(1.02/5)/Pmax
w=0.38
axC1.bar(units-w/2,np.full(5,0.852),w,color=BLUE)
bars=axC1.bar(units+w/2,r_mw,w,color='#d98880')
for k in np.where(r_mw>1)[0]: bars[k].set_color(RED)
axC1.axhline(1.0,color=RED,ls='--',lw=1.8)
axC1.text(4.45,1.05,'rating',color=RED,fontsize=9.5,ha='right')
axC1.annotate('pinned',(2.19,1.06),ha='center',fontsize=9.5,color=RED,fontweight='bold')
axC1.set_ylim(0,1.28); axC1.set_xticks(units); axC1.set_xticklabels([f'u{i+1}' for i in units],fontsize=9)
axC1.set_ylabel(r'$r_i$',fontsize=10.5); axC1.tick_params(labelsize=9)
axC1.spines[['top','right']].set_visible(False)
axC1.set_title('equal per-unit loading (blue)  vs  MW-averaging (red)',fontsize=9.8,color=DGRAY,pad=3)

axC2=fig.add_axes([0.665,0.150,0.315,0.295])
t=d['t']
axC2.plot(t,d['fF'],color=GREEN,lw=2.2)
axC2.plot(t,d['fI'],color=RED,lw=2.2)
axC2.axhline(49.0,color='#888888',ls='--',lw=1.4)
axC2.text(0.15,50.03,'feasible: recovers',fontsize=9.5,color=GREEN)
axC2.text(4.55,49.55,r'infeasible: $-0.233$ Hz/s',fontsize=9.5,color=RED)
axC2.text(7.9,49.06,'shed',fontsize=9,color='#777777',ha='right')
axC2.set_xlim(0,8); axC2.set_ylim(48.35,50.24)
axC2.set_xlabel('time (s)',fontsize=9.5,labelpad=4); axC2.set_ylabel('f (Hz)',fontsize=10)
axC2.tick_params(labelsize=9); axC2.spines[['top','right']].set_visible(False)

# ================= validation strip =================
fig.text(0.5,0.028,'Validated: 8 scenarios  ·  IEEE 14-bus  ·  delay, noise, packet loss, '
         r'$\pm$20% parameters, partition  ·  6 published methods on one identical testbed  ·  full code released',
         ha='center',fontsize=10.6,color='#333333',
         bbox=dict(boxstyle='round,pad=0.35',fc='#eef2f7',ec='#c8d4e2',lw=1.0))

fig.savefig('graphical_abstract.png',dpi=300,facecolor=BG)
fig.savefig('graphical_abstract.tiff',dpi=300,facecolor=BG,pil_kwargs={'compression':'tiff_lzw'})
fig.savefig('graphical_abstract.pdf',facecolor=BG)
print("v2 whole-story GA rendered")
