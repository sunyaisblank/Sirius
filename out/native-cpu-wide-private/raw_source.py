"""Derive only the pinned private raw16 sources; no compiler or workload mode."""
from pathlib import Path
import hashlib
import json
import re
import sys

R=Path(__file__).resolve().parents[2];W=Path(__file__).resolve().parent
B=R/'bin/linux-gcc/native-cpu-wide-private';S=B/'ispc-source'
ISPC=R/'bin/linux-gcc/native-cpu-simd-rows/compiler/bin/ispc'
TOOL_PIN=dict(bytes=211552768,sha256='1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e')
IDENTITY=R/'out/native-cpu-local-identities/identity_store.isph'
IDENTITY_PIN=dict(bytes=3052,sha256='0d464b6ce820abc38301d8f305ba6802ae04e5c13eff80c953bad118b786be2c')
PACKING=R/'out/native-cpu-inline-private-bank/packinghelper.isph'
PACKING_PIN=dict(bytes=1903,sha256='afce0e14e3f2732509c835e94d6178545dbd861d4f7cbd8dea59d462f5de173f')
CONTRACTS={
 'out/native-cpu-strict-ir-linking/canonical_private_dp_contract.json':dict(bytes=11905,sha256='242e1e101927473daea1c0cb8a9238d438e5f569a27d6218e01b56c6b46a674e'),
 'out/native-cpu-strict-ir-linking/canonical_private_remaining_contract.json':dict(bytes=58311,sha256='f2f9b0ec0fb0930485eaf0016d219327044f65602d22a9696b47e9450f98e1cc')}
PROFILES={
 'transport':dict(source='bin/linux-gcc/native-cpu-canonical-closed-transport/ispc-source/transport.ispc',bytes=780343,
  sha256='e3e55ac5f9a9cf8d32f3414ae29c4681702e07099cffd574a761980543c8e8b7',raw='transport',all_calls=3,tuples=3144,roots=40,layers=159),
 'endpoint':dict(source='bin/linux-gcc/native-cpu-canonical-closed-endpoint/ispc-source/endpoint.ispc',bytes=1186962,
  sha256='afafadf5a989ae1da57d2409e1855ef6852978596cdf3d8872e11e96e90e4956',raw='endpoint',all_calls=3,tuples=4148,roots=100,layers=195),
 'dopri_phase':dict(source='bin/linux-gcc/native-cpu-canonical-closed-dp/ispc-source/dopri_phase.ispc',bytes=710648,
  sha256='7e3dc7d0c602d7720248b5c8536e0d7b8ac295dc611afc77227b606f85cadef6',raw='dp',all_calls=4,tuples=2088,roots=200,layers=45)}
MATH=('dense_add','dense_multiply','dense_subtract','dense_divide','dense_sqrt',
      'tagged_add','tagged_multiply','tagged_subtract','tagged_divide','tagged_sqrt')


def pin(path):
    with path.open('rb') as stream:return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())
def identity(data):return dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
def doc(path):return json.loads(path.read_text())
def exact(path,expected):assert pin(path)=={k:expected[k] for k in ('bytes','sha256')},path


def apply(text,before,after,count,edits):
    assert text.count(before)==count,(before,text.count(before),count)
    assert after not in text,after
    edits.append(dict(before=before,after=after,count=count))
    return text.replace(before,after)


def inverse(text,edits,original):
    for edit in reversed(edits):
        assert text.count(edit['after'])==edit['count'],edit
        text=text.replace(edit['after'],edit['before'])
    assert text.encode()==original


def source_authority():
    for name,expected in CONTRACTS.items():exact(R/name,expected)
    dp=doc(R/next(name for name in CONTRACTS if 'private_dp_' in name))
    remaining=doc(R/next(name for name in CONTRACTS if 'remaining_' in name))
    assert dp['source_proof']['canonical_tuple_records_match_exact_admitted_table']
    assert dp['source_proof']['all_validation_before_first_publication']
    assert dp['source_proof']['sole_supported_caller_required200_lane0']
    assert dp['source_proof']['original_tuple_records']==2088
    assert dp['derived_source']=={k:PROFILES['dopri_phase'][k] for k in ('bytes','sha256')}
    for stage in ('transport','endpoint'):
        proof=remaining['source_proofs'][stage]
        for key in ('source_inverse_to_publication_exact','source_inverse_to_AS_exact',
            'canonical_tuple_records_match_exact_tables','mixed_families_supported',
            'all_root_validation_before_any_publication','original_combined_mask_five_word_store_order'):
            assert proof[key] is True,(stage,key)
        assert proof['original_tuple_records']==PROFILES[stage]['tuples']
        assert proof['original_layer_call_sites']==PROFILES[stage]['layers']
        assert remaining['derived_sources'][stage]=={k:PROFILES[stage][k] for k in ('bytes','sha256')}
    ep=remaining['source_proofs']['endpoint']
    for key in ('live_Endpoint11241_guard_unchanged','Endpoint_prefix_state_and_whole_state_rollback_unchanged',
                'required20_and100_supported_per_lane'):assert ep[key] is True,key
    exact(ISPC,TOOL_PIN);exact(IDENTITY,IDENTITY_PIN);exact(PACKING,PACKING_PIN)
    assert 'programCount' in PACKING.read_text()
    assert sum(p['tuples'] for p in PROFILES.values())==9380
    return dp,remaining


def planned():
    dp,remaining=source_authority();outputs={};proofs={};inputs={}
    paths=[Path(__file__).resolve(),W/'subgroup8.isph',IDENTITY,PACKING,ISPC,*[R/name for name in CONTRACTS],
        R/'out/native-cpu-value-fact-stores/value_facts.isph',R/'out/native-cpu-normalization-tags/normalization_tag.isph',
        R/'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html']
    for stage in ('transport','endpoint','dopri_phase'):
        paths.append(R/('bin/linux-gcc/native-cpu-slang-bridge/typed-batch/'+stage+'.program.hpp'))
    # Identity changes only its reduction domains and return variability.
    original=IDENTITY.read_bytes();text=original.decode();edits=[]
    text=apply(text,'static inline uniform bool SiriusTryIdentity_0(',
               'static inline varying bool SiriusTryIdentity_0(',1,edits)
    assert len(re.findall(r'\ball\(',text))==1 and len(re.findall(r'\bany\(',text))==3
    text=apply(text,'all(','SiriusSubgroup8All_0(',1,edits)
    text=apply(text,'any(','SiriusSubgroup8Any_0(',3,edits)
    inverse(text,edits,original)
    identity_path=S/'identity_subgroup8.isph';outputs[identity_path]=text.encode()
    identity_proof=dict(input=dict(path=str(IDENTITY.relative_to(R)),**IDENTITY_PIN),
        output=dict(path=str(identity_path.relative_to(R)),**identity(text.encode())),
        exact_edits=edits,source_inverse_exact=True,return_is_varying=True,old_all_calls=1,old_any_calls=3)
    for stage,profile in PROFILES.items():
        path=R/profile['source'];exact(path,profile);paths.append(path)
        original=path.read_bytes();text=original.decode();edits=[]
        old_export='sirius_simd_closed_'+profile['raw']+'_packet';new_export='sirius_wide_'+profile['raw']+'_packet'
        text=apply(text,'#define sirius_simd_'+stage+'_packet '+old_export,
                   '#define sirius_simd_'+stage+'_packet '+new_export,1,edits)
        for suffix in MATH:
            before='#define sirius_native_'+suffix+'8 sirius_native_v4_'+suffix+'8'
            after='#define sirius_native_'+suffix+'8 sirius_wide_pointer_'+suffix+'16'
            text=apply(text,before,after,1,edits)
        old_include='#include "'+str(IDENTITY)+'"'
        new_include='#include "'+str(W/'subgroup8.isph')+'"\n#include "'+str(identity_path)+'"'
        text=apply(text,old_include,new_include,1,edits)
        text=apply(text,'uniform uint32_t SiriusPackedBank[120];','uniform uint32_t SiriusPackedBank[240];',5,edits)
        text=apply(text,'uniform uint32_t SiriusPackedTags[24];','uniform uint32_t SiriusPackedTags[48];',5,edits)
        assert len(re.findall(r'\ball\(',text))==profile['all_calls'] and not re.search(r'\bany\(',text)
        text=apply(text,'all(','SiriusSubgroup8All_0(',profile['all_calls'],edits)
        inverse(text,edits,original)
        # The inverse is complete source equality; these semantic witnesses
        # independently retain the original literal operator/layer sequence.
        tuples=lambda s:re.findall(r'\bEvaluateOneStatic_0<\s*\d+U\s*,\s*\d+U\s*>\(\s*\d+U\s*,\s*\d+U\s*,\s*\d+U\s*,\s*\d+U\s*,\s*[A-Za-z_]\w*\s*,\s*[A-Za-z_]\w*\s*,\s*[A-Za-z_]\w*\s*\)',s)
        before_tuples,after_tuples=tuples(original.decode()),tuples(text)
        assert len(before_tuples)==len(after_tuples)==profile['tuples']
        assert after_tuples==before_tuples
        tuple_witness=identity('\n'.join(after_tuples).encode())
        for token in ('programCount','programIndex','for (uniform uint32_t invocation = 0; invocation < 64U; ++invocation)'):
            assert text.count(token)==original.decode().count(token),token
        assert text.count('*(&OUTPUT_WORD(')==original.decode().count('*(&OUTPUT_WORD(')
        target=S/(stage+'.ispc');outputs[target]=text.encode()
        proofs[stage]=dict(input=dict(path=profile['source'],**{k:profile[k] for k in ('bytes','sha256')}),
            output=dict(path=str(target.relative_to(R)),**identity(text.encode())),raw_export=new_export,
            exact_edits=edits,source_inverse_exact=True,original_tuple_count=profile['tuples'],
            literal_template_call_witness=dict(count=len(after_tuples),**tuple_witness),
            original_layer_call_sites=profile['layers'],root_positions=profile['roots'],
            canonical_tables_and_literal_operator_order_unchanged=True,all_validation_before_publication_unchanged=True,
            source_RHS_trials_quartics_and_refusals_unchanged=True,public_original_sources_untouched=True,
            programCount_indexing_and_full_context_initialization_unchanged=True,
            packed_payload_words=240,packed_tag_words=48,compiler_private_packing_bytes=1152,
            replaced_raw_all_calls=profile['all_calls'],replaced_raw_any_calls=0,
            reduction_scope='Three arithmetic shortcuts retain separate incoming low/high8 decisions; DP additionally preserves its canonical guard per subgroup.')
    for path in paths:inputs[str(path.relative_to(R))]=pin(path)
    commands={stage:[str(ISPC),str(S/(stage+'.ispc')),'-O2','--target=avx512skx-x16',
        '--cpu=skx','--arch=x86-64','--addressing=32','--PIC','--opt=disable-fma','--math-lib=system',
        '--emit-llvm-text','-o',str(B/(stage+'.ll'))] for stage in PROFILES}
    record=dict(source_only=True,compiler_run=False,source_inputs=inputs,
        outputs={str(path.relative_to(R)):identity(data) for path,data in outputs.items()},
        identity_header=identity_proof,stages=proofs,prospective_emit_commands=commands,
        programCount=16,old_math_core_lanes=8,canonical_tuple_total=9380,
        public_IO_ABI='Same six scalar pointer/count/first_x/groups arguments. New private symbols only; original public sources, buffers, Vulkan workgroup64 and logical IDs remain unchanged.',
        subgroup_proof='Each reducer captures incoming lanemask. Inactive predicate values are selected to the neutral value without freezing. Under that same mask, all(value|other-half) or any(value&same-half) produces a uniform result for each old8 subgroup; varying select returns that subgroup result. Identity return becomes varying so a completed/refused half retires alone; its remaining bytes and all data/tag snapshots/stores inverse exactly.',
        math_boundary='Ten new private pointer-helper names are prospective lowering inputs only. They must all be removed/bound by the reviewed strict/value16 lowering before native link/run. Existing math.o and all C++ arithmetic stay unchanged; root owns the9 value16 wrappers splitting into original8-lane cores. This source derivation claims no ABI lowering or zero-mask call proof beyond the explicit per8 raw reductions.',
        required_future_gates='Fresh pinned ordinary IR/intrinsic vocabulary, strict FP/mask/init/value16->math8 joins, context/staging/metadata bounds, async entry unwind/fullGP, v4 ISA/noFMA, private16/8/tail task partition and admission, current genuine-job C8/all32 and every original benefit/control/frame gate remain required. No candidate numerical/performance/qualification transfer.',
        scope='Private raw source-only width16 candidate. No compiler, math rewrite, Mesa/provider/scheduler mutation, workload, benchmark, extra warmup or scientific/release acceptance.')
    return outputs,record


def main():
    assert sys.argv[1:]==['derive'],'Only explicit source derivation is supported.'
    assert not S.exists() and not (W/'raw-source-contract.json').exists()
    outputs,record=planned()
    assert not any(path.exists() for path in outputs)
    S.mkdir(parents=True)
    for path,data in outputs.items():
        with path.open('xb') as stream:stream.write(data)
        assert pin(path)==record['outputs'][str(path.relative_to(R))]
    target=W/'raw-source-contract.json'
    with target.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(source_only=True,compiler_run=False,outputs=record['outputs'],contract=dict(path=str(target.relative_to(R)),**pin(target)))))


if __name__=='__main__':main()
