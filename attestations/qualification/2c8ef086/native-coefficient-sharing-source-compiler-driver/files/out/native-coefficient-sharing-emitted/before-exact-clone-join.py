"""Passive exact saved-SPIRV comparison; standard library only, no project imports.
Module normalization is preserved from the fd530 strict reader. No pure-CSE or
RTConstant exemption is imported. ROOT owns actual compiler and SPIR-V tools.
"""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import hashlib,json,re,difflib
R=Path(__file__).resolve().parents[2]
W=R/'out/native-coefficient-sharing-compiler'
P=Path(__file__).resolve().parent
BASE=R/'attestations/qualification/5b2df58b/native-constant-controls-benefit-pass/files/out/native-constant-compiler-source/candidate'
TOK=re.compile(r'"(?:\\.|[^"\\])*"|[^\s]+')
NAMES=['retained_transport']
def identity(p):
 b=Path(p).read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def pin(p):return {'path':str(p),**identity(p)}
def signature(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))

class Module:
    def __init__(self,path):
        self.path=path; self.text=path.read_text()
        self.defs=dict(re.findall(r'^\s*(%\S+) = (.+)$',self.text,re.M))
        self.names=dict(re.findall(r'OpName (%\S+) "([^"]*)"',self.text))
        self.bodies=dict(re.findall(r'^\s*(%\S+) = (OpFunction\b.*?OpFunctionEnd)',self.text,re.M|re.S))
        self.keys={};self.by_key={};self.rows={};self.kept={};self.local={};counts=Counter()
        for ident,body in self.bodies.items():
            name=self.names[ident];key=(name,counts[name]);counts[name]+=1
            self.keys[ident]=key;self.by_key[key]=ident
            rows=[ident+' = '+body.splitlines()[0]]+body.splitlines()[1:]
            self.rows[key]=[TOK.findall(s) for s in rows if s.strip()]
            selected=list(enumerate(self.rows[key]))
            self.kept[key]=selected
            for i,row in selected:
                if len(row)>1 and row[1]=='=':self.local[row[0]]=(key,len([v for v in self.local.values() if v[0]==key]))
        self.globals={i:e for i,e in self.defs.items() if e.startswith('OpVariable ') and e.split()[2]!='Function'}
        for i in self.globals:
            if i not in self.names:
                builtin=re.findall(r'OpDecorate '+re.escape(i)+r' BuiltIn (\S+)',self.text)
                assert len(builtin)==1 and builtin[0] in ('WorkgroupId','LocalInvocationIndex')
                self.names[i]='BuiltIn:'+builtin[0]
        assert len({self.names[i] for i in self.globals})==len(self.globals)
        self.static={i:e for i,e in self.defs.items() if i not in self.keys and i not in self.local and e.startswith(('OpType','OpConstant','OpSpecConstant','OpExtInstImport ','OpUndef '))}
        self.norm=lru_cache(None)(self._norm)

    def _norm(self,i):
        if i in self.keys:return ('function',self.keys[i])
        if i in self.globals:return ('global',self.names[i])
        if i in self.local:return ('local',self.local[i])
        assert i in self.static,(i,self.defs.get(i))
        return tuple(self.norm(v) if v.startswith('%') else v for v in TOK.findall(self.static[i]))

    def normalized_rows(self,key):
        return [tuple(self.norm(t) if t.startswith('%') else t for t in row) for _,row in self.kept[key]]

    def globals_normal(self):
        return {self.names[i]:tuple(self.norm(t) if t.startswith('%') else t for t in TOK.findall(e)) for i,e in self.globals.items()}

    def semantic_header(self):
        values=[]
        for line in self.text.splitlines():
            row=TOK.findall(line.strip())
            if not row:continue
            if row[0]=='OpEntryPoint':
                values.append(tuple(row[:2])+tuple([self.norm(row[2]),row[3],tuple(sorted(self.norm(t) for t in row[4:]))]))
            elif row[0] in ['OpCapability','OpExtension','OpMemoryModel','OpExecutionMode','OpExecutionModeId','OpDecorate','OpMemberDecorate','OpDecorateId']:
                values.append(tuple(self.norm(t) if t.startswith('%') else t for t in row))
        return Counter(values)


def main():
 findings=[]
 for name in NAMES:
  baseline=BASE/(name+'.final.spvasm');candidate=W/'candidate'/(name+'.final.spvasm')
  a=Module(baseline);b=Module(candidate)
  assert set(a.by_key)==set(b.by_key),('functions',name)
  assert a.globals_normal()==b.globals_normal(),('globals',name)
  diffs=[]
  for key in sorted(a.by_key):
   if key[0]=='CacheCoefficient':continue
   x,y=a.normalized_rows(key),b.normalized_rows(key)
   if x!=y:diffs.append({'function':list(key),'lengths':[len(x),len(y)],'first':next((i for i,(u,v) in enumerate(zip(x,y)) if u!=v),None)})
  extra=Counter(b.norm(i) for i in b.static)-Counter(a.norm(i) for i in a.static)
  missing=Counter(a.norm(i) for i in a.static)-Counter(b.norm(i) for i in b.static)
  findings.append({'name':name,'outside_differences':diffs,'extra_statics':[(k,v) for k,v in extra.items()],'missing_statics':[(k,v) for k,v in missing.items()],'helper_keys':[k for k in a.by_key if k[0]=='CacheCoefficient']})
 path=P/'strict-outside-finding.json';assert not path.exists();path.write_text(json.dumps(findings,indent=2)+'\n');print(json.dumps(findings))

if __name__=='__main__':main()
