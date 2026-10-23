"""Extend the frozen preservation adapter with actual frame disposition and exact reconstruction."""
from pathlib import Path
import hashlib
_adapter = Path(__file__).resolve().with_name('preserve.py')
assert hashlib.sha256(_adapter.read_bytes()).hexdigest() == 'ab88e7d21781bf63d2b27e767b8c37e3c0ac739172cbc7730f2e1e2e338097af'
_adapter_text = _adapter.read_text()
_tail = "exec(compile(_source,str(_basis),'exec'))\n"
assert _adapter_text.endswith(_tail)
exec(compile(_adapter_text[:-len(_tail)], str(_adapter), 'exec'))
_old = "frame.update(attempted=True,accepted=o['passed'],XML_present=(p.parent/'gtest.xml').exists())"
_new = """frame.update(attempted=True,owner_passed=o['passed'],accepted=False,XML_present=(p.parent/'gtest.xml').exists())
    disposition_path=W/'cold-frame/result.json'
    disposition=doc(disposition_path)
    assert disposition['attempted'] and disposition['source']==SOURCE
    assert disposition['stop_reason']==o['stop_reason'] and disposition['child_exit']==o['child_exit']
    assert disposition['wall_seconds']==o['wall_seconds']
    if disposition['accepted']:
        assert disposition['audit_completed'] and not disposition['failures'] and o['passed']
    frame.update(accepted=disposition['accepted'],disposition_identity=pin(disposition_path),disposition=disposition)"""
assert _source.count(_old) == 1
_source = _source.replace(_old, _new)
_old = 'Authenticate prior bundle; reconstruct authentic Mesa25.2.8 + exact5 IEEE32 prerequisite + original intergroup prototype, then apply the additive machine-codegen patch with zero fuzz. Source contract and ROOT admission bind all8 changed files and11100-file source identity. Configured original build inputs are sealed in the build owner manifests; actual staged provider is a content-addressed product.'
_new = 'Authenticate prior bundle; reconstruct authentic Mesa25.2.8 + exact5 IEEE32 prerequisite + original263 intergroup prototype, then apply bulk-zero.patch with zero fuzz. Source contract and ROOT admission bind the4 changed files and11100-file source identity. The rejected CodeGenNone patch is absent; original Default machine optimization is restored. Configured original build inputs are sealed in the build owner manifests; actual staged provider is a content-addressed product.'
assert _source.count(_old) == 1
_source = _source.replace(_old, _new)
exec(compile(_source, str(_basis), 'exec'))
