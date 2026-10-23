"""Complete the strict saved-code audit allowing only bijective local SSA renaming."""
from pathlib import Path
import hashlib
_basis=Path(__file__).resolve().with_name('read_code.py')
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='8d60d058a5b2b2f048b7db77d505d5758436880868f3656b6301ea9477e4a53c'
_wrapper=_basis.read_text()
_tail="exec(compile(_prefix+'\\n'+_setup+'\\n'+_loop+'\\n'+_end,str(_basis),'exec'))\n"
assert _wrapper.endswith(_tail)
exec(compile(_wrapper[:-len(_tail)],str(_basis),'exec'))
_addon=r'''
_strict_normalized_clear=normalized_clear
_alpha_rows=[]
def normalized_clear(ir,count,candidate):
    text=_strict_normalized_clear(ir,count,candidate)
    assert '%"' not in text, 'Quoted local identifiers require separate review'
    pattern=r'%[-A-Za-z$._0-9]+'
    definitions=re.findall(r'^  (%[-A-Za-z$._0-9]+) = ',text,re.M)
    assert len(definitions)==len(set(definitions))==62908
    names=set(re.findall(pattern,text))
    mapping={name:f'%sirius_alpha_{i}' for i,name in enumerate(definitions)}
    assert not set(mapping.values()) & names
    canonical=re.sub(pattern,lambda m:mapping.get(m[0],m[0]),text)
    _alpha_rows.append({'candidate':candidate,'definitions':definitions,
       'unmapped_names':sorted(names-set(definitions)),
       'canonical_sha256':hashlib.sha256(canonical.encode()).hexdigest()})
    return canonical
def alpha_proof():
    candidate,prior=_alpha_rows
    assert candidate['candidate'] is True and prior['candidate'] is False
    assert candidate['unmapped_names']==prior['unmapped_names']
    assert candidate['canonical_sha256']==prior['canonical_sha256']
    changed=[{'prior':a,'candidate':b} for a,b in zip(prior['definitions'],candidate['definitions']) if a!=b]
    assert len(changed)==193
    return {'definition_count_each':62908,'changed_definition_names':193,
      'changed_bijection':changed,'unmapped_arguments_types_and_block_names_unchanged':True,
      'remaining_IR_canonical_sha256':candidate['canonical_sha256'],
      'scope':'Global injective definition-order alpha renaming only. Every use replaced consistently and complete remaining text exact, including flags/constants/attributes/control flow. No Boolean or arithmetic rewrite; source/IR reasoning, not a machine-checked proof or numerical admission.'}
'''
_prefix+='\n'+_addon
_old='Remaining identifiers, operations, constants, flags, types and control flow exact.'
_new='Global bijective local SSA renaming permitted; every operation, constant, flag, type, attribute and control-flow instruction remains exact after substitution.'
assert _loop.count(_old)==1;_loop=_loop.replace(_old,_new)
_old="     'actual_provider':"
assert _end.count(_old)==1;_end=_end.replace(_old,"     'remaining_IR_alpha_proof':alpha_proof(),\n"+_old)
exec(compile(_prefix+'\n'+_setup+'\n'+_loop+'\n'+_end,str(_basis),'exec'))
