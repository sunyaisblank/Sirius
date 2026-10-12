"""Compare new observational metadata with immutable original producers/tables."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, struct

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASELINE = WORK/'baseline_retained_program.py'
CURRENT = ROOT/'src/sirius/kernels/retained_program.py'
HEADER = ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
MAPPING = ROOT/'attestations/review/9fbd0e8/stage7-write-mapping/files/mapping-result.json'
def seal(p):
    data=p.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def load(p, name):
    namespace={'__name__':name};exec(compile(p.read_bytes(),str(p),'exec'),namespace);return namespace
before={str(p.relative_to(ROOT)):seal(p) for p in (BASELINE,CURRENT,HEADER,MAPPING,Path(__file__))}
assert seal(BASELINE)['sha256']=='3754f0911b2b1a3d58276b23f4468d61a6092c65a1b522d7b557e373c6e96516'
assert seal(HEADER)['sha256']=='3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
assert seal(MAPPING)['sha256']=='6649a641dba12b2295476dd1951a33b7af190cd5ea2ff4a023c10cd40f8e6544'
old,new=load(BASELINE,'original_producer'),load(CURRENT,'candidate_producer')
arrays={}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',HEADER.read_text(),re.S):
    values=tuple(int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body))
    assert len(values)==int(count) and name not in arrays;arrays[name]=values
assert len(arrays)==40
def encode(kind,p):
    head=[p['instructions'],p['registers']]
    if kind not in ('Camera','RayCamera'):head.append(len(p['outputs']))
    return tuple(head+p['outputs']+p['operations']+p['layer_offsets'])
builders=[('Camera','build_camera_program'),('Transport','build_transport_program'),
          ('Endpoint','build_endpoint_program'),('Dense','build_dense_program'),
          ('Initialize','build_initialize_program'),('RayCamera','build_ray_camera_program'),
          ('DopriPhase','build_dopri_phase_program'),
          ('Transport','build_schwarzschild_transport_program'),
          ('Endpoint','build_schwarzschild_endpoint_program')]
defaults={};legacy={};matched=[]
for kind,name in builders:
    for parallel in (False,True):
        expected=old[name](parallel);actual=new[name](parallel)
        assert actual==expected,(name,parallel)
        assert 'source_plan' not in actual
        (defaults if parallel else legacy)[name]=actual
        matched.append({'builder':name,'parallel':parallel,'instructions':actual['instructions'],'registers':actual['registers']})
for kind,name in builders[:7]:
    encoded=encode(kind,defaults[name])
    if kind in ('Transport','Endpoint'):
        encoded+=encode(kind,defaults['build_schwarzschild_'+kind.lower()+'_program'])
    assert encoded==arrays['k'+kind+'Program'],kind
owner=new['build_stage_seven_owner_plan']()
reference=json.loads(MAPPING.read_text())
assert reference['accepted_source_mapping'] is True
assert reference['source_revision']=='9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
assert reference['runtime_equality_proved'] is False and reference['performance_trial_authorized'] is False
families={}
for family,old_family in (('general','general'),('schwarzschild','Schwarzschild')):
    data=owner[family]
    for stage in ('transport','endpoint'):
        name='build_'+('schwarzschild_' if family=='schwarzschild' else '')+stage+'_program'
        assert data['programs'][stage]==defaults[name]
        observed=new[name](True,with_plan=True)
        metadata=observed.pop('source_plan');assert observed==defaults[name]
        assert len(metadata['writes'])==observed['instructions']
    records=[r for r in reference['families'][old_family]['matches'] if r['opcode'] not in (0,1)]
    expected=tuple((r['Transport']['instruction_ordinal'],r['Endpoint']['instruction_ordinal']) for r in records)
    assert data['writes']==expected,family
    assert len(expected)==reference['families'][old_family]['common_arithmetic_nodes']
    assert len(expected)=={'general':71,'schwarzschild':26}[family]
    families[family]={'mapped_values':len(expected),'all_original_write_pairs_exact':True,
                      'prefix_values':sum(e<data['programs']['endpoint']['prefix_instructions'] for _,e in expected),
                      'suffix_values':sum(e>=data['programs']['endpoint']['prefix_instructions'] for _,e in expected)}
for name in ('build_transport_program','build_schwarzschild_transport_program','build_endpoint_program','build_schwarzschild_endpoint_program'):
    try:new[name](False,with_plan=True)
    except ValueError:pass
    else:raise AssertionError(('nonparallel metadata unexpectedly accepted',name))
assert before=={str(p.relative_to(ROOT)):seal(p) for p in (BASELINE,CURRENT,HEADER,MAPPING,Path(__file__))}
result={'accepted':True,'inputs':before,'default_builder_comparisons':matched,'original_program_arrays_exact':7,
        'original_header_arrays_retained':40,'families':families,'scope':'Observational producer metadata/default serialized program and original write correspondence only. No shader compilation, replay, numerical/performance or qualification verdict.'}
p=WORK/'plan-validation.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
