"""Join actual first field ABI/FP/metadata captures to the sealed canonical source."""
from pathlib import Path
import json
import re
import subprocess


def check(root,action,observation):
    slugs=('transport','endpoint','dopri_phase')
    symbols={}
    for line in subprocess.check_output(['/usr/bin/nm','-S','-a',action['runtime_provider']],text=True).splitlines():
        words=line.split()
        if len(words)==4 and words[2] in ('t','T'):
            symbols[words[3]]=int(words[0],16)
    for category,suffix,prefix in (('calls','_rows_try_packet','sirius_native_'),
        ('preflight_calls','_admit_packet','sirius_native_'),('field_calls','_field_try_packet','sirius_native_'),
        ('layer_calls','_field_layer','sirius_simd_')):
        rows=observation[category]
        assert len(rows)==3 and {row['function'] for row in rows}=={prefix+slug+suffix for slug in slugs}
    capability,=observation['capability_calls']
    assert capability['function']=='sirius_native_simd_supported' and capability['handled']
    for category in ('calls','preflight_calls','field_calls','layer_calls','capability_calls'):
        for row in observation[category]:
            assert row['returned'] and row['handled'] and row['thread']==row['return_thread']
            assert row['entry_sequence']<row['return_sequence'] and row['library']==action['runtime_provider']
            bases={int(line.split()[0].split('-')[0],16)-int(line.split()[2],16) for line in row['provider_maps']}
            assert row['pc']-symbols[row['function']] in bases
            if category in ('preflight_calls','capability_calls'):
                assert row['mxcsr']==row['return_mxcsr']
            else:
                assert row['mxcsr']&0xe040==row['return_mxcsr']&0xe040==0
            if category in ('calls','field_calls'):
                assert row['groups']==1 and row['first_x']<row['grid_x']<=row['capacity']<=65536
                assert row['y']==row['z']==0 and row['grid_y']==row['grid_z']==1
    schedule=json.loads((root/'out/native-cpu-field-simd/schedule-preparation.json').read_text())
    base=json.loads((root/'out/native-cpu-static-program/header-preparation.json').read_text())
    for row in observation['layer_calls']:
        slug=row['function'][len('sirius_simd_'):-len('_field_layer')]
        stage=base['headers'][slug]
        stride=stage['input_row_words_from_original_packet_guard']
        output_stride={'transport':4699,'endpoint':3829,'dopri_phase':4149}[slug]
        words=[int(v) for v in re.findall(r'(\d+)u',(root/stage['table_path']).read_text())]
        assert row['input_words']==1+row['capacity']*stride+len(words)
        assert row['output_words']==row['capacity']*output_stride
        assert row['row']<row['capacity'] and row['row_base']==row['row']*output_stride
        family,=[item for item in schedule['families'][slug]['programs'] if row['program_base']==1+row['capacity']*stride+item['base']]
        assert row['count']==family['instructions'] and row['registers']==family['registers'] and row['layer']<family['layers']
        offset=family['base']+3+family['roots']+5*family['instructions']+row['layer']
        assert row['start']==words[offset] and row['end']==words[offset+1]
        assert row['scratch'] and row['status'] and row['status_before']==row['status_after']==0
        assert any(item['function'] and 'sirius_generated_'+slug+'_field::' in item['function'] for item in row['callers'])
    namespace={};path=root/'out/native-cpu-field-cohort/compatibility.py'
    exec(compile(path.read_text(),str(path),'exec'),namespace)
    component=namespace['check'](root)['scope']
    return {'first_rows_returns':3,'first_pure_returns':3,'first_field_control_returns':3,
        'first_field_layer_returns':3,'runtime_ISA_supported':True,'canonical_layer_metadata_and_private_context_captured':True,
        'transport_component':component,
        'scope':'Actual first entry/return per category/stage, provider address/mapping, canonical layer args and original FP scope. Independent concurrent first calls are not claimed as one shared packet. All original32complete readbacks remain separate mandatory checks.'}
