"""Read the original two finite probe results without repeating dispatch."""
from pathlib import Path
import hashlib,json,xml.etree.ElementTree as ET
W=Path(__file__).resolve().parent
oracle=json.loads((W/'upward-oracle.json').read_text())
expected=(W/'upward-expected.bin').read_bytes(); variants={}
for label in ['baseline','candidate']:
    folder=W/('primitive-'+label); owner=json.loads((folder/'owner.json').read_text())
    records=[json.loads(line) for line in (folder/'stdout.log').read_text().splitlines()]
    assert owner['passed'] and owner['provider_mappings'] and owner['whole_inputs_unchanged']
    selected=records[0]; actual=records[-1]
    assert selected['device']=='llvmpipe (LLVM 20.1.2, 256 bits)'
    assert selected['kind']=='software' and selected['driver']=='llvmpipe'
    assert selected['driver_info']=='Mesa 25.2.8-0ubuntu0.24.04.2 (LLVM 20.1.2)'
    assert selected['RTE32']==1
    # dladdr uses /lib while proc maps and the sealed ICD use /usr/lib on this
    # merged-/usr installation. Resolve the actual path, preserving byte seals.
    assert Path(selected['resident_driver_path']).resolve(strict=True)==Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so').resolve(strict=True)
    assert actual['completed_dispatches']==1 and actual['cases']==oracle['cases']
    assert actual['observed_words']==oracle['words_per_variant']
    assert actual['mismatches']==0 and actual['untouched_complement_words']==0
    assert actual['actual_resident_bytes']<=8*1024*1024
    assert (W/('primitive-'+label+'-actual.bin')).read_bytes()==expected
    variants[label]={'selected':selected,'actual':actual,'owner_wall_seconds':owner['wall_seconds']}
xml=ET.parse(W/'admission/gtest.xml').getroot()
assert xml.get('tests')=='3'
assert all(int(xml.get(key,'0'))==0 for key in ['failures','errors','skipped','disabled'])
assert len(xml.findall('.//testcase'))==3 and not any(xml.findall('.//'+tag) for tag in ['failure','error','skipped'])
for case in xml.findall('.//testcase'):
    assert case.get('status')=='run' and case.get('result')=='completed'
result={'pass':True,'source_revision':json.loads((W/'retention-gate-frozen.json').read_text())['candidate_revision'],
    'cases_per_variant':oracle['cases'],'words_per_variant':oracle['words_per_variant'],'classes':oracle['classes'],
    'halfway_parities':oracle['admitted_halfway_parities'],'variants':variants,
    'expected_sha256':hashlib.sha256(expected).hexdigest(),
    'scope':'Finite actual RPUpMultiply Vulkan invocation against independent exact-Fraction words, old and candidate; selected software provider and full output canaries verified. Not universal proof, primitive performance or frame/science qualification.'}
(W/'primitive-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'pass':True,'cases_per_variant':oracle['cases'],'words_per_variant':oracle['words_per_variant'],
    'mismatches':0,'original_admission_tests':3}))
