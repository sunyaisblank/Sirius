"""Independent defining metric, matrix inversion and complete-family differentiation.
Finite high precision references are witnesses, not certified intervals.
"""
from pathlib import Path
import struct,json
import mpmath as mp

def f32(x):return struct.unpack('<f',struct.pack('<f',x))[0]
def fixtures():
 names=['flat_pinhole','flat_moving_pupil','schwarzschild_pupil','schwarzschild_moving_pinhole','kerr_positive_pupil','kerr_negative_pupil','charged_kerr_pupil','cosmological_pupil']
 families=[(0,0,0,0),(0,0,0,0),(1,0,0,0),(1,0,0,0),(1,.7,0,0),(1,-.7,0,0),(1,.5,.2,0),(1,0,0,.001)]
 rows=[]
 for i,family in enumerate(families):
  thin=i not in [0,3];beta=[0,0,0] if i in [0,2] else [.125,.0625,-.03125]
  row=[*family,.25,4,1,.5,-1,.05,.02,.01,0,-1,0,1,.03,*beta,1920,1080,f32(1920/1080),.625,10,701.25,413.5,.125 if thin else 0,-.0625 if thin else 0,int(thin),i,0]
  assert len(row)==32;rows.append([f32(x) for x in row])
 return names,rows

def metric(row,x):
 M,a,Q,L=row[:4];g=mp.diag([-1,1,1,1])
 if M==a==Q==L==0:return g
 xx,y,z=x[1:];reduced=xx*xx+y*y+z*z-a*a
 r=mp.sqrt((reduced+mp.sqrt(reduced*reduced+4*a*a*z*z))/2)
 H=(2*M*r-Q*Q)/(r*r+a*a*z*z/(r*r))
 if a==0:H+=L*r*r/3
 ell=mp.matrix([1,(r*xx+a*y)/(r*r+a*a),(r*y-a*xx)/(r*r+a*a),z/r])
 return g+H*ell*ell.T

def dot(g,a,b):return (a.T*g*b)[0]
def frame(row,x):
 g=metric(row,x);inv=g**-1;lapse=1/mp.sqrt(-inv[0,0]);u=mp.matrix([-lapse*inv[i,0] for i in range(4)])
 axes=[mp.matrix([0,*row[8:11]]),mp.matrix([0,*[-v for v in row[11:14]]]),mp.matrix([0,*row[14:17]])];basis=[]
 for raw in axes:
  e=raw+u*dot(g,raw,u)
  for old in basis:e-=old*dot(g,e,old)
  basis.append(e/mp.sqrt(dot(g,e,e)))
 beta=row[17:20];b2=sum(v*v for v in beta)
 if b2==0:return [u,*basis]
 gamma=1/mp.sqrt(1-b2);bv=sum((basis[i]*beta[i] for i in range(3)),mp.zeros(4,1));coef=(gamma-1)/b2
 return [gamma*(u+bv),*[basis[i]+coef*beta[i]*bv+gamma*beta[i]*u for i in range(3)]]

def family(row,z):
 central=mp.matrix(row[4:8]);initial=frame(row,central);fx,fy,pr,pu=z;thin=row[29]==1
 x=central+initial[3]*pr+initial[2]*pu if thin else central
 F=frame(row,x);D=row[24] if thin else mp.mpf(1);T=row[23]
 q=mp.matrix([D,D*T*(1-2*fy/row[21])-(pu if thin else 0),D*T*(2*fx/row[20]-1)*row[22]-(pr if thin else 0)])
 n=q/mp.sqrt(sum(v*v for v in q));k=-F[0]+sum((F[i+1]*n[i] for i in range(3)),mp.zeros(4,1))
 return list(x)+list(k)+[v for axis in F for v in axis]

def calculate(raw,dps):
 with mp.workdps(dps):
  row=list(map(mp.mpf,raw));z=row[25:29];base=family(row,z);x=mp.matrix(base[:4]);k=mp.matrix(base[4:8]);u=mp.matrix(base[8:12]);g=metric(row,x);inv=g**-1
  def derivative_at(function,args,c):
   return mp.diff(lambda t:function(args[:c]+[t]+args[c+1:]),args[c])
  deriv=[[derivative_at(lambda zz:family(row,zz)[j],z,c) for j in range(12)] for c in range(4)]
  dg=[mp.matrix(4,4) for _ in range(4)]
  for c in range(1,4):
   for i in range(4):
    for j in range(4):dg[c][i,j]=derivative_at(lambda xx:metric(row,xx)[i,j],list(x),c)
  G=[[[sum(inv[m,s]*(dg[a][s,b]+dg[b][s,a]-dg[s][a,b]) for s in range(4))/2 for b in range(4)] for a in range(4)] for m in range(4)]
  X=[r[:4] for r in deriv];K=[r[4:8] for r in deriv];du=[r[8:12] for r in deriv]
  V=[[K[c][m]+sum(G[m][a][b]*k[a]*X[c][b] for a in range(4) for b in range(4)) for m in range(4)] for c in range(4)]
  null=[];freq=[]
  for c in range(4):
   v=mp.matrix(V[c]);nabla=mp.matrix([du[c][m]+sum(G[m][a][b]*u[a]*X[c][b] for a in range(4) for b in range(4)) for m in range(4)])
   null.append(dot(g,k,v));freq.append(dot(g,v,u)+dot(g,k,nabla))
  scientific=base+[v for group in [X,K,V,du] for column in group for v in column]+base[8:24]
  assert len(scientific)==104
  return {'scientific':[mp.nstr(v,dps) for v in scientific],'null':[mp.nstr(v,dps) for v in null],'frequency':[mp.nstr(v,dps) for v in freq],'metric':[[mp.nstr(g[i,j],dps) for j in range(4)] for i in range(4)],'Gamma':[[[mp.nstr(v,dps) for v in a] for a in m] for m in G]}


def prepare():
    output = Path(__file__).resolve().parents[3] / "out" / "retained-camera-build"
    output.mkdir(parents=True, exist_ok=True)
    names, rows = fixtures()
    references = []
    for name, row in zip(names, rows):
        low, high = calculate(row, 100), calculate(row, 180)
        with mp.workdps(190):
            stability = max(abs(mp.mpf(a) - mp.mpf(b)) / (1 + abs(mp.mpf(b)))
                            for a, b in zip(low["scientific"], high["scientific"]))
            assert stability < mp.mpf("1e-90")
            assert max(abs(mp.mpf(value)) for key in ["null", "frequency"]
                       for value in high[key]) < mp.mpf("1e-160")
        references.append({"name": name, "input": row, "reference100": low,
                           "reference180": high, "stability": str(stability)})
        print(name + ": independent reference identities passed", flush=True)
    (output / "moderate-camera-inputs.bin").write_bytes(
        b"".join(struct.pack("<32f", *row) for row in rows))
    (output / "moderate-camera-references.json").write_text(
        json.dumps(references, indent=2) + "\n")


if __name__ == "__main__":
    prepare()
