"""Read-only full-entry CFI/actual GP-save machine proof, object and final DSO."""
from pathlib import Path
import bisect
import hashlib
import importlib.util
import re
import subprocess
import sys
R=Path.cwd();sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'))
from elf_proof import Elf
GP=dict(rbx='c-56',rbp='c-16',r12='c-48',r13='c-40',r14='c-32',r15='c-24')
PROFILES={'transport':dict(offset=0,stack='40ec0300',returns=[1371]),
          'endpoint':dict(offset=27,stack='c0680300',returns=[1503])}


def verify(path,symbol,source_table=None):
    path=Path(path)
    if symbol=='sirius_closed_dp_packet':
        authority=R/'out/native-cpu-closed-dp-unwind/verify.py'
        assert hashlib.sha256(authority.read_bytes()).hexdigest()=='b4bb26d693e8a4c500011e3eb7ca2e8fa70b389eb3a0d32960b3a24a6f3c1d00'
        spec=importlib.util.spec_from_file_location('AX_exact_DP_unwind',authority)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.Elf=Elf
        result=module.verify_unwind(path,symbol)
    else:
        stage,=[s for s in PROFILES if symbol=='sirius_closed_'+s+'_packet'];profile=PROFILES[stage]
        elf=Elf(path);function,=[s for s in elf.symbols if s['name']==symbol and s['section'] and s['info']&15==2 and s['size']]
        start,end=function['value'],function['value']+function['size'];code=elf.function_bytes(function);p=profile['offset']
        prologue=bytes.fromhex('554889e54157415641554154534883e4c04881ec'+profile['stack'])
        assert code[p:p+24]==prologue
        assert code[:64].count(bytes.fromhex('554889e5415741564155415453'))==1
        # These push/mov/and/sub bytes prove the real frame and six saved GP
        # slots; the two stage profiles come from their actual native products.
        frames=subprocess.check_output(['/usr/bin/readelf','-W','--debug-dump=frames-interp',str(path)],text=True)
        fdes=list(re.finditer(r'^([0-9a-f]+) .*\bFDE .*\bpc=([0-9a-f]+)\.\.([0-9a-f]+)\s*$',frames,re.M))
        i,=[i for i,m in enumerate(fdes) if (int(m[2],16),int(m[3],16))==(start,end)]
        block=frames[fdes[i].end():fdes[i+1].start() if i+1<len(fdes) else len(frames)]
        columns=None;rows=[]
        for line in block.splitlines():
            if line.strip().startswith('LOC '):assert columns is None;columns=line.split()
            elif columns and re.match(r'^[0-9a-f]+\s',line):
                fields=line.split();assert len(fields)==len(columns);row=dict(zip(columns,fields))
                offset=int(row.pop('LOC'),16)-start;assert 0<=offset<=function['size'];rows.append(dict(offset=offset,**row))
            elif rows and not line.strip():break
        assert rows and rows[0]['offset']==0 and all(a['offset']<b['offset'] for a,b in zip(rows,rows[1:]))
        offsets=[r['offset'] for r in rows]
        def state(offset):
            j=bisect.bisect_right(offsets,offset)-1;assert j>=0;return rows[j]
        assert state(0)['CFA']==state(p)['CFA']=='rsp+8' and state(0)['ra']=='c-8'
        assert state(p+1)['CFA']=='rsp+16' and state(p+1)['rbp']=='c-16'
        assert state(p+4)['CFA']=='rbp+16'
        assert state(p+24)['CFA']=='rbp+16' and all(state(p+24).get(reg)==rule for reg,rule in GP.items())
        assembly=subprocess.check_output(['/usr/bin/objdump','-d','--no-show-raw-insn','--start-address='+str(start),'--stop-address='+str(end),str(path)],text=True)
        instructions=[]
        for line in assembly.splitlines():
            m=re.match(r'^\s*([0-9a-f]+):\s+([a-z][a-z0-9.]*)\b(?:\s+(.*))?$',line)
            if m:instructions.append((int(m[1],16)-start,m[2],m[3] or ''))
        calls=[off for off,op,args in instructions if op in ('call','callq')]
        returns=[off for off,op,args in instructions if op in ('ret','retq')]
        assert calls and returns==profile['returns']
        for off in calls:
            rules=state(off);assert rules['CFA']=='rbp+16' and rules['ra']=='c-8'
            assert all(rules.get(reg)==rule for reg,rule in GP.items()),(symbol,off,rules)
        epilogue=bytes.fromhex('488d65d85b415c415d415e415f5dc5f877c3')
        for off in returns:
            assert code[off+1-len(epilogue):off+1]==epilogue
            assert state(off)['CFA']=='rsp+8' and state(off)['ra']=='c-8'
        assert all(row['ra']=='c-8' and row['CFA'] in ('rsp+8','rsp+16','rbp+16') for row in rows)
        result=dict(full_function_FDE=True,symbol=symbol,function_bytes=function['size'],
            relative_range=[0,function['size']],relative_CFI_rows=rows,initial_CFA='rsp+8',body_CFA='rbp+16',
            complete_GP_save_rules=True,saved_GP_rules=GP,saved_return_PC='CFA-8',shrinkwrap_and_epilogues_checked=True,
            actual_GP_push_and_restore_machine_bytes_checked=True,prologue_offset=p,prologue_machine_bytes=prologue.hex(),
            epilogue_machine_bytes=epilogue.hex(),covered_call_sites=len(calls),checked_return_offsets=returns,
            same_thread_live_job_and_return_observation_still_required=True,
            scope='Exact current stage prologue/GP push slots and restored epilogue bytes, full FDE/CFA/GP rules at every call and return. No DP-layout assumption or runtime unwind/numerical verdict.')
    if source_table:
        for key in ('relative_range','relative_CFI_rows','covered_call_sites','checked_return_offsets'):
            assert result[key]==source_table[key],(symbol,key)
    return result
