import json,importlib.util
from pathlib import Path
import mpmath as mp
mp.mp.dps=75
spec=importlib.util.spec_from_file_location('camera_reference',str(Path(__file__).resolve().parent.parent / 'retained_camera/reference.py'));ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
reflection=[-1,1,-1,1];params=list(map(mp.mpf,[1,.998,0,0]))
def metric(x):
 g=ref.metric(params,mp.matrix([x[i]*reflection[i] for i in range(4)]))
 return mp.matrix([[g[i,j]*reflection[i]*reflection[j] for j in range(4)] for i in range(4)])
def geometry(x,second=True):
 g=metric(x);inv=g**-1;dg=[mp.zeros(4) for _ in range(4)];dd=[[mp.zeros(4) for _ in range(4)] for _ in range(4)]
 for a in range(4):
  for b in range(a,4):
   f=lambda *v:metric([x[0],*v])[a,b]
   for axis in range(1,4):
    n=[0,0,0];n[axis-1]=1;dg[axis][a,b]=dg[axis][b,a]=mp.diff(f,tuple(x[1:]),tuple(n))
    if second:
     for other in range(axis,4):
      order=n[:];order[other-1]+=1
      value=mp.diff(f,tuple(x[1:]),tuple(order));dd[axis][other][a,b]=dd[axis][other][b,a]=dd[other][axis][a,b]=dd[other][axis][b,a]=value
 G=[[[sum(inv[m,s]*(dg[a][s,b]+dg[b][s,a]-dg[s][a,b]) for s in range(4))/2 for b in range(4)] for a in range(4)] for m in range(4)]
 return g,inv,dg,dd,G
def mat(v):return mp.matrix(list(map(mp.mpf,v)))
def pack(x,p,columns):return mp.matrix(list(x)+list(p)+[v for X,P in columns for v in [*X,*P]])
def unpack(y):return y[:4,0],y[4:8,0],[(y[8+8*c:12+8*c,0],y[12+8*c:16+8*c,0]) for c in range(4)]
def initialize(row):
 x=mat(row['x']);k=mat(row['k']);g,inv,dg,dd,G=geometry(list(x),False);columns=[]
 for c in row['columns']:
  X,V=mat(c['X']),mat(c['V']);K=mp.matrix([V[m]-sum(G[m][a][b]*k[a]*X[b] for a in range(4) for b in range(4)) for m in range(4)])
  P=g*K+mp.matrix([sum(dg[a][m,b]*X[a]*k[b] for a in range(4) for b in range(4)) for m in range(4)])
  columns.append((X,P))
 return pack(x,g*k,columns)
def rhs(y):
 x,p,columns=unpack(y);g,inv,dg,dd,G=geometry(list(x));k=inv*p
 force=mp.matrix([(k.T*dg[m]*k)[0]/2 for m in range(4)])
 result=[]
 for X,P in columns:
  K=inv*(P-sum((dg[a]*k*X[a] for a in range(4)),mp.zeros(4,1)))
  dP=mp.matrix([(K.T*dg[m]*k)[0]+sum((k.T*dd[a][m]*k)[0]*X[a]/2 for a in range(4)) for m in range(4)])
  result.append((K,dP))
 return pack(k,force,result)
def rk4(y,h):
 a=rhs(y);b=rhs(y+a*h/2);c=rhs(y+b*h/2);d=rhs(y+c*h)
 return y+(a+2*b+2*c+d)*h/6

def project(y):
 x,p,columns=unpack(y);g,inv,dg,dd,G=geometry(list(x),False);k=inv*p
 # Differentiate the temporal quadratic implicitly, keeping all spatial
 # coordinate tangent derivatives. This is independent of the covariant
 # first-kind reduction used by the C++ experiment.
 A=g[0,0];B=2*sum(g[0,i]*k[i] for i in range(1,4));C=sum(g[i,j]*k[i]*k[j] for i in range(1,4) for j in range(1,4))
 roots=[(-B+sgn*mp.sqrt(B*B-4*A*C))/(2*A) for sgn in [-1,1]]
 kp=k.copy();kp[0]=min(roots,key=lambda v:abs(v-k[0]));out=[]
 for X,P in columns:
  K=inv*(P-sum((dg[a]*k*X[a] for a in range(4)),mp.zeros(4,1)))
  numerator=sum((kp.T*dg[a]*kp)[0]*X[a] for a in range(4))+2*sum(g[a,b]*kp[a]*K[b] for a in range(4) for b in range(1,4))
  K[0]=-numerator/(2*sum(g[0,a]*kp[a] for a in range(4)))
  V=K+mp.matrix([sum(G[m][a][b]*kp[a]*X[b] for a in range(4) for b in range(4)) for m in range(4)])
  out.append({'X':[str(v) for v in X],'V':[str(v) for v in V]})
 return {'x':[str(v) for v in x],'k':[str(v) for v in kp],'columns':out}

def freeze_metric_cases():
    rows=json.loads(Path(__file__).with_name('inputs.json').read_text())['states'];out=[]
    for i,row in enumerate(rows):
     samples=[]
     for precision in [75,110]:
      mp.mp.dps=precision
      g,inv,dg,dd,G=geometry(list(map(mp.mpf,row['previous']['x'])),False)
      samples.append([g[a,b] for a in range(4) for b in range(4)]+[inv[a,b] for a in range(4) for b in range(4)]+[dg[axis][a,b] for axis in range(4) for a in range(4) for b in range(4)])
     fields=[]
     for low,high in zip(*samples):
      hi=float(high);lo=float(high-mp.mpf(hi));fields.append({'hi':hi.hex(),'lo':lo.hex(),'reference':str(high),'precision_gap':str(abs(low-high))})
     out.append({'position':row['previous']['x'],'fields':fields});print(i,flush=True)
    Path(__file__).with_name('metric_reference.json').write_text(json.dumps({'precision':[75,110],'description':'Independent defining Kerr metric, generic inverse and numerical high-precision differentiation; observed precision gap is not an interval certificate.','cases':out},indent=2)+'\n')
    lines=['// Generated independent precision witnesses; regenerate with reference.py.', '// clang-format off', '#pragma once','#include <array>','namespace sirius::test::critical_fixture {','struct Pair {double hi,lo;};','struct MetricCase {std::array<double,4> position;std::array<Pair,96> fields;};','inline constexpr std::array<MetricCase,10> metric_cases{{']
    for r in out:
     lines.append('    {{'+','.join(float(v).hex() for v in r['position'])+'}, {{')
     lines.extend('        {'+v['hi']+','+v['lo']+'},' for v in r['fields'])
     lines.append('    }}},')
    lines+=['}};','} // namespace sirius::test::critical_fixture','// clang-format on','']
    Path(__file__).with_name('metric_reference.h').write_text('\n'.join(lines))

def freeze_transport_cases():
    mp.mp.dps=75
    folder=Path(__file__).resolve().parent
    rows=json.loads((folder/'inputs.json').read_text())['states']
    inputs=[];references=[]
    for index,row in enumerate(rows):
        projected=project(initialize(row['previous']))
        initial={key:list(map(float,projected[key])) for key in ['x','k']}
        initial['columns']=[{field:list(map(float,column[field])) for field in ['X','V']} for column in projected['columns']]
        inputs.append({'previous':initial,'h':row['h']})
        state=initialize(initial);h=mp.mpf(row['h'])
        full=project(rk4(state,h));refined=project(rk4(rk4(state,h/2),h/2))
        gap=max(abs(mp.mpf(full['columns'][column][field][axis])-mp.mpf(refined['columns'][column][field][axis])) for column in range(4) for field in ['X','V'] for axis in range(4))
        references.append({'index':index,'gap':str(gap),'full':full,'refined':refined})
        print(index,'RK4 refinement gap',gap,flush=True)
    (folder/'projected_inputs.json').write_text(json.dumps({'description':'Derived one-step oracle fixtures. The ten original readbacks are independently projected onto the null constraint at 75 digits, including the derivative, then materialized once as binary64 inputs. This does not replace or alter original_launches.json or inputs.json.','states':inputs},indent=2)+'\n')
    (folder/'transport_reference.json').write_text(json.dumps({'description':'Independent 75-digit RK4 one-step and two-half-step refinement at each derived input. Observed refinement gap, not an interval certificate.','cases':references},indent=2)+'\n')
    def values(row):return row['x']+row['k']+[value for column in row['columns'] for field in ['X','V'] for value in column[field]]
    lines=['// Independent one-step transport witnesses. See README.md.','// clang-format off','#pragma once','#include <array>','namespace sirius::test::critical_fixture {','struct TransportCase {double h;std::array<double,40> initial;std::array<double,40> expected;double refinement_gap;};','inline constexpr std::array<TransportCase,10> transport_cases{{']
    for initial,reference in zip(inputs,references):
        lines.append('    {'+float(initial['h']).hex()+', {')
        lines.extend('        '+float(value).hex()+',' for value in values(initial['previous']))
        lines.append('    }, {')
        lines.extend('        '+float(value).hex()+',' for value in values(reference['refined']))
        lines.append('    }, '+float(reference['gap']).hex()+'},')
    lines+=['}};','} // namespace sirius::test::critical_fixture','// clang-format on','']
    (folder/'transport_reference.h').write_text('\n'.join(lines))

def freeze_launch_header():
    folder=Path(__file__).resolve().parent
    launches=json.loads((folder/'original_launches.json').read_text())['launches']
    assert len(launches)==14
    lines=['// Exact saved launch values; do not regenerate camera samples to replace them.','// clang-format off','#pragma once','#include <array>','namespace sirius::test::critical_fixture {','struct Launch {int x,y;double u,v;std::array<double,4> direction;std::array<double,3> boost;};','inline constexpr std::array<Launch,14> original_launches{{']
    for row in launches:
        assert len(row)==6 and len(row[4])==4 and len(row[5])==3
        lines.append('    {'+str(row[0])+','+str(row[1])+','+float(row[2]).hex()+','+float(row[3]).hex()+', {'+','.join(float(value).hex() for value in row[4])+'}, {'+','.join(float(value).hex() for value in row[5])+'}},')
    lines+=['}};','} // namespace sirius::test::critical_fixture','// clang-format on','']
    (folder/'original_launches.h').write_text('\n'.join(lines))

if __name__ == '__main__':
 import argparse
 parser=argparse.ArgumentParser(description='Explicit independent critical-transport fixture regeneration (mpmath required).')
 parser.add_argument('--metric-fixtures', action='store_true')
 parser.add_argument('--transport-fixtures', action='store_true')
 parser.add_argument('--launch-header', action='store_true')
 args=parser.parse_args()
 if args.metric_fixtures:freeze_metric_cases()
 if args.transport_fixtures:freeze_transport_cases()
 if args.launch_header:freeze_launch_header()
 if not args.metric_fixtures and not args.transport_fixtures and not args.launch_header:parser.print_help()
