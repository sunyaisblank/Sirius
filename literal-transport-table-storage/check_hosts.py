"""Recheck complete private ELF spans and their read-only load segments."""
from pathlib import Path
import hashlib, json, struct
WORK = Path(__file__).resolve().parent

def main():
    preparation = json.loads((WORK / 'host-preparation.json').read_text())
    assert preparation['pass'] and preparation['all_inputs_unchanged']
    owner = json.loads((WORK / 'prepare-hosts/owner.json').read_text())
    assert owner['passed'] and owner['remaining_owned_processes'] == [] and not owner['cleanup_errors']
    inputs = json.loads((WORK / 'host-inputs-before.json').read_text())
    for name, item in inputs.items():
        p = Path(name); raw = p.read_bytes()
        assert len(raw) == item['bytes'] and hashlib.sha256(raw).hexdigest() == item['sha256'], name
    report = {'pass': True, 'inputs_rechecked': len(inputs), 'scope': 'ELF symbol-address byte spans, allocated non-write/non-execute sections and readable non-write/non-execute PT_LOAD segments; no runtime mapping or execution conclusion', 'hosts': {}}
    for mode, artifact in preparation['artifacts'].items():
        exe = Path(artifact['executable']['path']); raw = exe.read_bytes()
        assert len(raw) == artifact['executable']['bytes'] and hashlib.sha256(raw).hexdigest() == artifact['executable']['sha256']
        assert raw[:6] == b'\x7fELF\x02\x01'
        phoff, shoff = struct.unpack_from('<QQ', raw, 32)
        phsize, phnum, shsize, shnum = struct.unpack_from('<HHHH', raw, 54)
        assert phsize == 56 and shsize == 64 and phoff + phsize*phnum <= len(raw) and shoff + shsize*shnum <= len(raw)
        programs = [struct.unpack_from('<IIQQQQQQ', raw, phoff + phsize*i) for i in range(phnum)]
        sections = [struct.unpack_from('<IIQQQQIIQQ', raw, shoff + shsize*i) for i in range(shnum)]
        records = {}
        for group in ('all40_private_spans', 'original36_test_reference_spans'):
            for name, item in artifact[group].items():
                at, address, size = item['file_offset'], item['address'], item['bytes']
                assert hashlib.sha256(raw[at:at+size]).hexdigest() == item['sha256']
                sections_here = [s for s in sections if s[1] != 8 and s[4] <= at and at+size <= s[4]+s[5] and s[3]+at-s[4] == address]
                assert len(sections_here) == 1 and sections_here[0][2] & 2 and not sections_here[0][2] & 5, (mode,name)
                loads = [p for p in programs if p[0] == 1 and p[2] <= at and at+size <= p[2]+p[5] and p[3]+at-p[2] == address]
                assert len(loads) == 1 and loads[0][1] == 4, (mode,name)
                records[group + '/' + name] = {'bytes': size, 'sha256': item['sha256'], 'section_flags': sections_here[0][2], 'load_flags': loads[0][1]}
        assert len(records) == 76
        report['hosts'][mode] = {'executable': artifact['executable'], 'spans': records}
    (WORK / 'host-payload-check.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'pass': True, 'inputs_rechecked': len(inputs), 'complete_readonly_spans': 152}))

if __name__ == '__main__':
    main()
