"""Compiler-only full-table canonical read substitution; no product or device write."""
from pathlib import Path
import hashlib, importlib.util, json, re, shutil, subprocess, sys
sys.dont_write_bytecode = True
ROOT=Path(__file__).resolve().parents[2]
WORK=Path(__file__).resolve().parent
SOURCE=ROOT/'src/sirius/kernels'
STAGE=WORK/'kernels'

def seal(path):
    raw=path.read_bytes()
    return {'path':str(path.relative_to(ROOT)), 'bytes':len(raw), 'sha256':hashlib.sha256(raw).hexdigest()}

def write(name, value):
    (WORK/name).write_text(json.dumps(value,indent=2)+'\n')

def words(program):
    return [program['instructions'],program['registers'],len(program['outputs']),*program['outputs'],*program['operations'],*program['layer_offsets']]

def replace(text, old, new, count=1):
    assert text.count(old)==count,(old,text.count(old),count)
    return text.replace(old,new)

def main():
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert revision=='b81e061fa081cdaed37eca3add950bd254bf296d'
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT)
    STAGE.mkdir()
    for p in sorted(SOURCE.glob('*.slang')): shutil.copyfile(p,STAGE/p.name)
    shutil.copyfile(SOURCE/'portable_binary32.h',STAGE/'portable_binary32.h')
    spec=importlib.util.spec_from_file_location('compact_transport_program',SOURCE/'retained_program.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    programs={'General':m.build_transport_program(parallel=True),'Schwarzschild':m.build_schwarzschild_transport_program(parallel=True)}
    encoded={n:words(p) for n,p in programs.items()}
    header=ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
    raw=re.search(r'kTransportProgram\{\{(.*?)\}\};',header.read_text(),re.S)[1]
    cached=[int(v,0) for v in re.findall(r'(0x[0-9a-fA-F]+|\d+)u?',raw)]
    assert encoded['General']+encoded['Schwarzschild']==cached and len(cached)==15967
    schedule=[]
    for family,p in programs.items():
        assert p['registers']==459 and len(p['outputs'])==40
        assert p['layer_offsets'][0]==0 and p['layer_offsets'][-1]==p['instructions']
        for layer,(start,end) in enumerate(zip(p['layer_offsets'],p['layer_offsets'][1:])):
            assert 0<end-start<=64
            dest,reads=set(),set()
            for index in range(start,end):
                op,d,a,b,c=p['operations'][5*index:5*index+5]
                assert 0<=op<=11 and 0<=d<459 and d not in dest;dest.add(d)
                if op==1: assert a<45
                elif op>=2:
                    assert a<459;reads.add(a)
                    if op not in (6,7,8): assert b<459;reads.add(b)
                    if op in (10,11): assert c<459;reads.add(c)
            assert not dest.intersection(reads)
            for lane in range(min(16,end-start)):
                schedule.append({'family':family,'layer':layer,'lane':lane,'nodes':list(range(start+lane,end,16))})
        assert sorted(i for s in schedule if s['family']==family for i in s['nodes'])==list(range(p['instructions']))
    original=(SOURCE/'retained_transport.slang').read_text()
    generic=replace(original,'bool EvaluateProgram(','bool EvaluateProgramGeneric(')
    one=generic[generic.index('bool EvaluateOne('):generic.index('bool EvaluateProgramGeneric(')]
    program=generic[generic.index('bool EvaluateProgramGeneric('):generic.index('RetainedTriple Rational(')]
    # A reverse substitution must recover each original function byte for byte.
    edits_one=[('bool EvaluateOne(', 'bool EvaluateOneCanonical('),
               ('uint programBase','bool schwarzschild'),('programBase + 43','43'),
               ('inputs[address]','CanonicalTransportWord(address, schwarzschild)'),
               *[(f'inputs[address + {i}]',f'CanonicalTransportWord(address + {i}, schwarzschild)') for i in range(1,5)]]
    canonical_one=one
    for a,b in edits_one: canonical_one=replace(canonical_one,a,b)
    restored=canonical_one
    for a,b in reversed(edits_one): restored=replace(restored,b,a)
    assert restored==one
    edits_program=[('bool EvaluateProgramGeneric(', 'bool EvaluateProgramCanonical('),
                   ('uint programBase','bool schwarzschild'),('programBase + 43','43'),
                   ('inputs[layerBase]', 'CanonicalTransportWord(layerBase, schwarzschild)'),
                   ('inputs[layerBase + layers]', 'CanonicalTransportWord(layerBase + layers, schwarzschild)'),
                   ('inputs[layerBase + layer]', 'CanonicalTransportWord(layerBase + layer, schwarzschild)'),
                   ('inputs[layerBase + layer + 1]', 'CanonicalTransportWord(layerBase + layer + 1, schwarzschild)')]
    canonical_program=program
    for a,b in edits_program: canonical_program=replace(canonical_program,a,b)
    for a,b,count in [('EvaluateOne(', 'EvaluateOneCanonical(',3),
                      ('rowBase, programBase, registers','rowBase, schwarzschild, registers',3),
                      ('inputs[programBase + 3 + i]','CanonicalTransportWord(3 + i, schwarzschild)',4)]:
        canonical_program=replace(canonical_program,a,b,count);edits_program.append((a,b))
    restored=canonical_program
    for a,b in reversed(edits_program): restored=restored.replace(b,a)
    assert restored==program and 'programBase' not in canonical_program and 'inputs[' not in canonical_one+canonical_program
    generated=[]
    for family,values in encoded.items():
        generated += [f'static uint Literal{family}Table[{len(values)}] = {{',
                      *[', '.join(f'{v}u' for v in values[i:i+16])+',' for i in range(0,len(values),16)],'};']
    generated += ['groupshared uint literalTransportMatch;',
                  'bool RecognizeLiteralTransport(uint programBase, bool schwarzschild, uint lane) {',
                  'if (lane == 0) literalTransportMatch = 1;','GroupMemoryBarrierWithGroupSync();',
                  'if (schwarzschild) {','for (uint i = lane; i < 2871u; i += 16u)',
                  '    if (inputs[programBase + i] != LiteralSchwarzschildTable[i])',
                  '        InterlockedAnd(literalTransportMatch, 0u);','} else {',
                  'for (uint i = lane; i < 13096u; i += 16u)',
                  '    if (inputs[programBase + i] != LiteralGeneralTable[i])',
                  '        InterlockedAnd(literalTransportMatch, 0u);','}',
                  'GroupMemoryBarrierWithGroupSync();','return literalTransportMatch != 0;','}',
                  'uint CanonicalTransportWord(uint offset, bool schwarzschild) {',
                  'if (schwarzschild) return LiteralSchwarzschildTable[offset];',
                  'return LiteralGeneralTable[offset];','}',canonical_one,canonical_program]
    block='\n'.join(generated)+'\n'
    modified=replace(generic,'RetainedTriple Rational(',block+'RetainedTriple Rational(')
    anchor='#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)\n    // The immutable plan leaves never-assigned slots invalid across all stages.'
    recognize='bool literalSchwarzschild = programBase != 1u + 230u * rows;\n    bool literalProgram = RecognizeLiteralTransport(programBase, literalSchwarzschild, lane);\n'
    modified=replace(modified,anchor,recognize+anchor)
    before='if (!EvaluateProgram(row, rowBase, programBase, count, registers, layers, stage, lane)) return;'
    router='if (literalProgram) {\n            if (!EvaluateProgramCanonical(row, rowBase, literalSchwarzschild, count, registers, layers, stage, lane)) return;\n        } else {\n            if (!EvaluateProgramGeneric(row, rowBase, programBase, count, registers, layers, stage, lane)) return;\n        }'
    modified=replace(modified,before,router)
    guard='#if !defined(SIRIUS_RETAINED_PORTABLE) || !defined(SIRIUS_RETAINED_NORMAL_SUM32) || !defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT) || SIRIUS_RETAINED_EXECUTION_LANES != 16\n#error This compiler-only diagnostic requires the existing portable normal Transport layout.\n#endif\n'
    modified=guard+modified
    inverse=modified.removeprefix(guard).replace(block,'',1).replace(recognize,'',1).replace(router,before,1).replace('bool EvaluateProgramGeneric(','bool EvaluateProgram(',1)
    assert inverse==original and 'outputs[rowBase + 3]' not in modified
    (STAGE/'retained_transport.slang').write_text(modified)
    write('programs.json',programs);write('schedule.json',schedule)
    write('source-frame-check.json',{'pass':True,'original_source':seal(SOURCE/'retained_transport.slang'),'candidate_source':seal(STAGE/'retained_transport.slang'),'canonical_one_inverse_byte_exact':True,'canonical_program_inverse_byte_exact':True,'whole_source_inverse_byte_exact':True,'original_generic_body_unchanged_except_declaration_name':True,'node_word_substitutions':5,'layer_word_substitutions':4,'root_word_substitutions':4,'marker_absent':True,'scope':'Only canonical program-word reads change. Numerical/scientific/compiler/runtime equivalence not implied.'})
    facts={'scope':'compiler-only; no device, numerical, performance or product adoption','revision':revision,'header':seal(header),'inputs':[seal(p) for p in sorted(SOURCE.glob('*.slang'))]+[seal(SOURCE/'retained_program.py'),seal(ROOT/'scripts/build-retained-kernels.py'),seal(SOURCE/'portable_binary32.h')],
           'programs':{n:{'words':len(encoded[n]),'nodes':p['instructions'],'layers':len(p['layer_offsets'])-1,'roots':40,'registers':459} for n,p in programs.items()},'constructor_words_equal_current_embedded_table':True,'recognition':'Exact selected full table after unchanged original admission/flat/selection; uniform match selects canonical indexed evaluator, mismatch unchanged Generic; no diagnostic output marker','emission':seal(STAGE/'retained_transport.slang'),'schedule':seal(WORK/'schedule.json')}
    write('emission.json',facts);print(json.dumps({'programs':facts['programs'],'source':facts['emission']}))

if __name__=='__main__': main()
