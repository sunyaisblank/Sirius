"""Read actual native word captures; separate control acceptance from equality."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import collections,hashlib,json,re,struct,xml.etree.ElementTree as ET
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
assert len(sys.argv)==2 and sys.argv[1] in ('device-timestamps','mixed-transport','factored-endpoints')
case=sys.argv[1];folder=WORK/case
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    data=p.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def words(data):
    assert len(data)%4==0;return struct.unpack('<%dI'%(len(data)//4),data)
def props(p):
    root=ET.parse(p).getroot();cases=list(root.iter('testcase'))
    assert len(cases)==1 and not any(list(root.iter(n)) for n in ('failure','error','skipped'))
    nodes=list(cases[0].iter('property'));assert len({n.get('name') for n in nodes})==len(nodes)
    return cases[0],{n.get('name'):n.get('value') for n in nodes}
plan=read(WORK/'execution-plan.json');schema=read(WORK/'value-schema.json')
assert schema['version']==1 and schema['row_words']==768
assert read(folder/'native-control-bridge.json')['accepted']
assert read(folder/'retirement-run.json')['accepted'] and read(folder/'native-retirement.json')['pass']
owner=read(folder/'run/owner.json');assert owner['numerical_control_passed'] and owner['source_revision']==plan['source_revision']
test,properties=props(folder/'run/gtest.xml')
assert test.get('classname')+'.'+test.get('name')==plan['cases'][case]
prior=ROOT/'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files'/case/'run'
old_paths=sorted((prior/'shader-clock-observations').glob('*.json'),key=lambda p:read(p)['ordinal'])
capture_props={k:json.loads(v) for k,v in properties.items() if k.startswith('value_witness_capture_')}
directory=folder/'run/word-value-observations'
paths=sorted(directory.glob('*.json'),key=lambda p:read(p)['ordinal'])
assert len(paths)==len(capture_props)==len(old_paths)>0
assert len(list(directory.glob('*.bin')))==3*len(paths)
assert {p.stem for p in paths}=={p.stem for p in directory.glob('*.bin') if not p.name.endswith(('.input.bin','.output.bin'))}

header=WORK/'baseline-products/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert seal(header)['sha256']=='3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
arrays={}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (k(?:Transport|Endpoint)Program)\{\{(.*?)\}\};',header.read_text(),re.S):
    values=tuple(int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body));assert len(values)==int(count)
    arrays[name]=values
assert len(arrays)==2
lengths={'transport':{'general':13096,'Schwarzschild':2871},'endpoint':{'general':11240,'Schwarzschild':9903}}
prefix={'general':250,'Schwarzschild':202}
key_kinds={'transport':{'general':[5,5,5,5,2,6,6,6,5],'Schwarzschild':[5,1,1,1,2,6,6,6,5]},
           'endpoint':{'general':[5,5,5,5,1,5,5,5,5],'Schwarzschild':[5,1,1,1,1,5,5,5,5]}}
ordinals=set();captures=[];totals=collections.Counter()
for path,old_path in zip(paths,old_paths):
    meta=read(path);old=read(old_path);ordinal=meta['ordinal']
    assert ordinal not in ordinals and meta==capture_props['value_witness_capture_'+str(ordinal)];ordinals.add(ordinal)
    for key in ('ordinal','stage','capacity','active_rows','groups_y','groups_z','slot','paired'):
        assert meta[key]==old[key],('Original capture topology changed',key,meta,old)
    assert meta['version']==1 and meta['stage'] in ('transport','endpoint')
    stage=meta['stage'];title=stage.title();capacity=meta['capacity'];active=meta['active_rows']
    assert 0<active<=capacity<=32 and meta['groups_y']==meta['groups_z']==1
    assert meta['words_per_row']==768 and meta['slot'] in (0,1) and (meta['paired'] or meta['slot']==0)
    assert meta['sidecar'] not in (meta['input'],meta['output'])
    assert 0<meta['sidecar_allocation_bytes']<=8*1024*1024
    assert meta['schema_sha256']==seal(WORK/'value-schema.json')['sha256']
    module=plan['exact_observer_modules']['retained_'+stage+'.spv'];assert meta['module_sha256']==module['sha256']
    expected_input=230 if stage=='transport' else 225;expected_output=4699 if stage=='transport' else 3829
    assert meta['input_row_words']==expected_input and meta['output_row_words']==expected_output
    buffers={}
    for suffix,size_key,hash_key in [('.bin','bytes','sha256'),('.input.bin','input_bytes','input_sha256'),('.output.bin','output_bytes','output_sha256')]:
        p=path.with_suffix(suffix);data=p.read_bytes();assert len(data)==meta[size_key] and hashlib.sha256(data).hexdigest()==meta[hash_key]
        buffers[suffix]=words(data)
    raw=buffers['.bin'];inp=buffers['.input.bin'];out=buffers['.output.bin'];canonical=arrays['k'+title+'Program']
    assert len(raw)==capacity*768 and len(inp)==1+capacity*expected_input+len(canonical) and len(out)==capacity*expected_output
    assert inp[0]==meta['stored_capacity']
    assert not any(raw[active*768:]),'Inactive sidecar rows were not cleared'
    base=1+inp[0]*expected_input;layout=inp[0]==capacity and base+len(canonical)==len(inp)
    supplied=inp[base:] if layout else ()
    combined=layout and supplied==canonical
    table_identity={}
    for family,length in lengths[stage].items():
        offset=schema['families'][title][family]['selected_table_offset_words']
        table_identity[family]=layout and supplied[offset:offset+length]==canonical[offset:offset+length]
    rows=[]
    for row in range(active):
        record=raw[row*768:(row+1)*768];output=out[row*expected_output:(row+1)*expected_output]
        keys=[record[16+7*i:23+7*i] for i in range(9)]
        arithmetic=[record[96+9*i:105+9*i] for i in range(71)]
        assert not any(record[7:16]) and not any(record[79:96]) and not any(record[735:])
        for i,key in enumerate(keys):
            if key[0]:assert key[0]&~7==0 and key[1]==schema['key_records'][i]
            else:assert not any(key)
        for value in arithmetic:
            assert value[0] in (0,1)
            if not value[0]:assert not any(value)
        family=None
        if record[0]:
            assert record[0]==1 and record[1]==(1 if stage=='transport' else 2)
            assert record[2] in (1,2);family='general' if record[2]==1 else 'Schwarzschild'
            assert record[3]==base+schema['families'][title][family]['selected_table_offset_words']
            assert record[4]==inp[0]
            if stage=='transport':assert record[5]==40 and record[6]==1
            else:assert (record[5],record[6]) in ((20,2),(100,6))
        else:assert not any(record[:7])
        canonical_selected=family is not None and table_identity[family]
        final_success=output[0]==1
        complete=canonical_selected and final_success and ((stage=='transport' and output[1:3]==(7,1)) or (stage=='endpoint' and record[5:7]==(100,6)))
        if canonical_selected:
            mapping=schema['families'][title][family]['arithmetic'];assert len(mapping) in (71,26)
            for index,value in enumerate(arithmetic):
                if index>=len(mapping):assert not any(value);continue
                if not value[0]:continue
                m=mapping[index];kind=6 if stage=='transport' else (20 if m['ordinal']<prefix[family] else 100)
                assert value[1:4]==(m['ordinal'],m['destination'],(kind<<8)|m['opcode'])
            if complete:
                assert [k[0] for k in keys]==key_kinds[stage][family]
                assert all(v[0]==1 for v in arithmetic[:len(mapping)])
        item={'row':row,'header':list(record[:7]),'original_output_header':list(output[:4]),'family_from_publication':family,'canonical_selected_table':canonical_selected,'original_final_success':final_success,'complete_word_witness':complete,'keys':[list(k) for k in keys],'arithmetic':[list(v) for v in arithmetic],'input_packet':list(inp[1+row*expected_input:1+(row+1)*expected_input]),'fifth_phase':list(output[4:204]) if stage=='transport' else None}
        rows.append(item);totals[stage+'_active_rows']+=1;totals[stage+'_complete_rows']+=int(complete)
        totals[stage+'_unclassified_or_noncanonical_rows']+=int(not canonical_selected)
        totals[stage+'_reached_arithmetic']+=sum(v[0] for v in arithmetic)
        totals[stage+'_reached_invalid_arithmetic']+=sum(v[0]==1 and v[8]==0 for v in arithmetic)
    captures.append({'capture':meta,'combined_table_canonical':combined,'selected_table_identities':table_identity,'layout_available':layout,'rows':rows})
    totals[stage+'_captures']+=1

stream=None
if case=='device-timestamps':
    _,prior_props=props(prior/'gtest.xml')
    current={k:v for k,v in properties.items() if '_readback_' in k}
    assert len(current)==32 and current=={k:v for k,v in prior_props.items() if '_readback_' in k}
    modules={k:v for k,v in properties.items() if k.startswith('selected_module_')}
    assert len(modules)==7 and set(modules)=={k for k in prior_props if k.startswith('selected_module_')}
    for index,stage in [(1,'transport'),(2,'endpoint')]:
        m=plan['exact_observer_modules']['retained_'+stage+'.spv'];assert modules['selected_module_'+str(index)]=='stage='+stage+';bytes='+str(m['bytes'])+';sha256='+m['sha256']
    assert all(v==prior_props[k] for k,v in modules.items() if k not in ('selected_module_1','selected_module_2'))
    marker=lambda v:dict(part.split('=',1) for part in v.split(';'))
    marks={k:marker(v) for k,v in properties.items() if k.startswith('device_observation_')}
    old_marks={k:marker(v) for k,v in prior_props.items() if k.startswith('device_observation_')}
    assert len(marks)==51 and set(marks)==set(old_marks)
    for k,v in marks.items():
        assert set(v)==set(old_marks[k])
        assert all(v[n]==old_marks[k][n] for n in ('phase','stage','available0','available1'))
        assert v['available0']==v['available1']=='1'
    stream={'original_complete_buffer_digests':32,'ordered_available_markers':51,'selected_modules':7,'exact_changed_modules':['transport','endpoint'],'prior_xml':seal(prior/'gtest.xml')}

pairs=[];pair_totals=collections.Counter();pending=None
if case=='device-timestamps':
    for capture in captures:
        meta=capture['capture']
        if meta['stage']=='transport':pending=capture;continue
        previous=pending;pending=None
        if previous is None or meta['active_rows']!=2*previous['capture']['active_rows']:continue
        for t,e in zip(previous['rows'],capture['rows'][:previous['capture']['active_rows']]):
            packet=t['input_packet'][:20]+t['fifth_phase']+t['input_packet'][220:225]
            raw_mismatch=[i for i,(a,b) in enumerate(zip(packet,e['input_packet'])) if a!=b]
            key_unobserved=[schema['key_records'][i] for i,(a,b) in enumerate(zip(t['keys'],e['keys'])) if not a[0] or not b[0]]
            key_mismatch=[schema['key_records'][i] for i,(a,b) in enumerate(zip(t['keys'],e['keys'])) if a[0] and b[0] and a[2:7]!=b[2:7]]
            eligible=not raw_mismatch and t['complete_word_witness'] and e['complete_word_witness'] and t['family_from_publication']==e['family_from_publication']
            exact_keys=eligible and not key_unobserved and not key_mismatch
            family=t['family_from_publication'];compared=[]
            if eligible:
                for index,m in enumerate(schema['families']['Endpoint'][family]['arithmetic']):
                    a=t['arithmetic'][index];b=e['arithmetic'][index]
                    differing=[i for i,(x,y) in enumerate(zip(a[4:9],b[4:9])) if x!=y]
                    compared.append({'slot':index,'kind':b[3]>>8,'different_word_indices':differing,'transport_words':a[4:9],'endpoint_words':b[4:9]})
                    pair_totals['arithmetic_records_compared']+=1
                    pair_totals['prefix_records' if b[3]>>8==20 else 'suffix_records']+=1
                    pair_totals['arithmetic_mismatches']+=int(bool(differing))
                    pair_totals['exact_key_arithmetic_records_compared']+=int(exact_keys)
                    pair_totals['exact_key_arithmetic_mismatches']+=int(exact_keys and bool(differing))
                pair_totals[family+'_complete_high_pairs']+=1
            pair_totals['immediate_high_pairs']+=1
            pair_totals['raw_packet_mismatches']+=int(bool(raw_mismatch))
            pair_totals['normalized_key_mismatches']+=int(bool(key_mismatch))
            pair_totals['pairs_with_unobserved_keys']+=int(bool(key_unobserved))
            pair_totals['unpaired_complete_word_witness']+=int(not eligible)
            pairs.append({'transport_ordinal':previous['capture']['ordinal'],'endpoint_ordinal':meta['ordinal'],'row':t['row'],'family':family,'raw_packet_different_indices':raw_mismatch,'normalized_key_different_records':key_mismatch,'unobserved_key_records':key_unobserved,'complete_canonical_high_pair':eligible,'exact_normalized_key_condition':exact_keys,'arithmetic_comparisons':compared})

result={'accepted_control_and_capture_analysis':True,'source_revision':plan['source_revision'],'case':plan['cases'][case],'counts':dict(totals),'stream':stream,'captures':captures,'immediate_high_pairs':pairs,'pair_counts':dict(pair_totals),'scope':'Finite actual original native controls and observer word records only. Pairing consumes only immediately following high rows, preserving separate lower/interior/trial calls. Equality is reported independently of control acceptance; invalid/reached/unpublished/noncanonical rows are retained. No universal native-word theorem, safe reuse, complete geometry ownership, ordinary benefit, frame or full qualification.'}
destination=folder/'word-analysis.json';assert not destination.exists();destination.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('accepted_control_and_capture_analysis','source_revision','case','counts','pair_counts')},indent=2))
