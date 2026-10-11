"""Recover only the native counted prefix; leave the rejected bridge unchanged."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def seal(p):return {'bytes':p.stat().st_size,'sha256':sha(p)}
out=WORK/'console-readiness-requested';bridge=WORK/'console-readiness-requested-bridge.json'
b=read(bridge);o=read(out/'owner.json');e=read(WORK/'console-readiness-requested-emergency.json')
assert not b['accepted'] and b['returncode']==0 and b['errors']==['AssertionError()'] and not b['outer_stop']
assert b['linux_child_birth_absent'] and b['emergency_observer']['linux_birth_absent'] and b['emergency_observer']['returncode']==0
assert o['status']=='completed' and not o['outer_stop'] and not o['cleanup_errors']
for k in ('console_ready_observed','observation_only_intentional_tool_stop','console_detached','console_geometry_verified_after_attachment','owned_process_absent','tool_unchanged','owned_service_absent','service_unchanged','developer_census_completed','intentional_service_stop'):assert o[k] is True,k
assert not o['accessible_remaining_developer_processes'] and o['returncode']==o['service_returncode']==2**32-1
assert e['recorded_births_absent'] and not e['errors']
assert o['outer_guard_seconds']==30 and o['rss_guard_bytes']==1073741824
assert all(seal(ROOT/p)==v for p,v in b['input_seals'].items())
cs=(WORK/'owned_console_requested.cs').read_text()
assert 'read!=length)throw Error();' in cs and 'string raw=cells.ToString();' in cs
observations=[]
for name,suffix in (('initial',0xe220),('ready',ord('x'))):
 p=out/(name+'-screen.json');s=read(p);n=s['width']*s['rows_read'];raw=s['cells']
 assert s['width']==o['required_console_width']==184 and s['height']==2048 and s['rows_read']==s['cursor_y']+1
 assert 0<=s['cursor_x']<184 and 0<=s['cursor_y']<2047 and len(raw)==n+1 and ord(raw[n])==suffix
 assert set(s['attached_pids'])=={o['owner_pid'],o['test_pid']} and len(s['attached_pids'])==2
 prefix=raw[:n];assert '\0' not in prefix
 lines=[prefix[i:i+184] for i in range(0,n,184)];assert all(line[-1]==' ' for line in lines)
 assert s['text']==''.join(line.rstrip(' ')+'\r\n' for line in lines)
 if name=='ready':
  for m in ('Profiling (RGP) feature enabled (auto-capture enabled)','Waiting for GPU application to connect...'):assert m in s['text']
 observations.append({'snapshot':str(p.relative_to(ROOT)),'snapshot_seal':seal(p),'native_requested_count':n,'native_read_count_source_inferred':n,'excluded_serialized_suffix_codepoint':suffix,'bounded_prefix_sha256_utf16le':hashlib.sha256(prefix.encode('utf-16le')).hexdigest(),'geometry':[s['width'],s['height']],'cursor':[s['cursor_x'],s['cursor_y']]})
assert not (out/'unused.rgp').exists()
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==b['source_revision']==subprocess.check_output(['git','rev-parse','@{u}'],text=True).strip()
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
record={'bounded_readiness_prefix_verified':True,'original_bridge_accepted':False,'original_bridge_seal':seal(bridge),'source_revision':b['source_revision'],'observations':observations,'recorded_ownership_closed':True,'all_frozen_inputs_exact':len(b['input_seals']),'native_count_limitation':'Count inferred from the immutable executed C# branch requiring read==length; not separately serialized. Excluded suffixes are unaccounted marshaling artifacts, not console cells or NUL terminators.','scope':'Separate derived counted-prefix observation only. Original rejection unchanged. No science, dispatch, trace, executed cost, performance benefit, frame or full qualification.','independent_reviewer':'/root/command_owner_review','references':['https://learn.microsoft.com/en-us/windows/console/readconsoleoutputcharacter','https://learn.microsoft.com/en-us/dotnet/standard/native-interop/best-practices']}
p=WORK/'console-requested-prefix-result.json';assert not p.exists();p.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'bounded_readiness_prefix_verified':True,'all_frozen_inputs_exact':len(b['input_seals']),'snapshots':2}))
